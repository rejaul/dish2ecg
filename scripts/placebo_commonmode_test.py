"""Is the ECGRDVQ-vs-ECGCIPA shift carried by the placebo group?

ECGCIPA has ONE placebo group (10 subjects, part 1) and its mean change from
baseline is subtracted, time-matched, from every active arm - the four
parallel arms *and* the dofetilide crossover arm. Any sampling error in that
one group therefore moves all ECGCIPA drug effects the same way. ECGRDVQ
subtracts each subject's own placebo period.

For each shared drug and study this script reports
  * the mean placebo correction actually applied to the arm's records
    (``dQTcF - ddQTcF``, ms);
  * the drug-only effect at a matched concentration, with NO placebo
    subtraction (baseline-corrected change only).
If the placebo group carries the shift, the drug-only effects agree between
studies and the applied placebo corrections differ by about the shift.

It also repeats the ECGRDVQ J-Tpeakc analysis with the power-law correction,
because the ECGRDVQ/ECGDMMLD reference script uses the linear form.
"""
from __future__ import annotations
import argparse, sys, warnings
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parent))
import clinical_reference as cr
import ecgcipa_reference as er
from dish2ecg import (fit_concentration_response, fridericia, jtpeakc_power,
                      subtract_baseline, double_delta)
from dish2ecg.io import load_physionet_clinical, load_ecgcipa
from dish2ecg.provenance import Endpoint

SHARED = {"Ranolazine": "Ranolazine", "Verapamil": "Verapamil HCL",
          "Dofetilide": "Dofetilide"}


def _fit(frame, col, ccol, gcol, at, prov, baseline_only=False):
    ep = Endpoint(values=frame[col].astype(float), name=col, units="ms",
                  domain="clinical", provenance=prov)
    # baseline_only=True is a DELIBERATE diagnostic: no placebo subtraction, so the
    # guardrail that normally demands control_subtraction is relaxed here, visibly.
    req = ("rate_correction", "baseline_subtraction") if baseline_only else (
        "rate_correction", "baseline_subtraction", "control_subtraction")
    f = fit_concentration_response(frame, ep, concentration_col=ccol,
                                   group_col=gcol, at_concentration=at,
                                   require_steps=req)
    return f.effect_at_concentration, f.effect_ci90[0], f.effect_ci90[1]


def main(out: Path) -> None:
    warnings.filterwarnings("ignore")
    out.mkdir(parents=True, exist_ok=True)

    # ---- ECGRDVQ (paired placebo, crossover) --------------------------------
    raw = load_physionet_clinical("ecgrdvq")
    d, _ = cr.prepare(raw, "ecgrdvq")
    g = cr.collapse_replicates(d).reset_index(drop=True)
    eps, dd, _, _ = cr.corrected_endpoints(g)
    stepR = subtract_baseline(g, eps["QTcF"], subject_col="RANDID", arm_col="arm",
                              baseline_mask=g["is_baseline"], match_on_time=False)
    g["dQTcF"] = stepR.values.to_numpy()
    g["ddQTcF"] = dd["QTcF"].values.to_numpy()
    g["plc_applied"] = g["dQTcF"] - g["ddQTcF"]

    # ---- ECGCIPA (one unpaired placebo group) -------------------------------
    D = load_ecgcipa()
    a = er.analysis_set(D["adeg"])
    qtcf = fridericia(a["QT"], a["RR"])
    stepC = subtract_baseline(a, qtcf, subject_col="USUBJID", arm_col="TRTA",
                              baseline_mask=a["is_baseline"], match_on_time=False)
    a["dQTcF"] = stepC.values.to_numpy()
    from dish2ecg import subtract_control
    ddC = subtract_control(a, stepC, arm_col="TRTA", control_arm="Placebo",
                           time_col="ATPTN", paired=False, name="ddQTcF")
    a["ddQTcF"] = ddC.values.to_numpy()
    a["plc_applied"] = a["dQTcF"] - a["ddQTcF"]

    rows = []
    for cipa_arm, rdvq_arm in SHARED.items():
        pk = er.concentrations_ngml(D["adpc"], er.ARM_ANALYTE[cipa_arm])
        s = a[(a["TRTA"] == cipa_arm) & (~a["is_baseline"])].merge(
            pk, on=["USUBJID", "ATPTN"], how="inner")
        s = s[s["conc"].notna() & s["ddQTcF"].notna() & s["dQTcF"].notna()]
        o = g[(g["arm"] == rdvq_arm) & g["conc_ngml"].notna() & (g["TPT"] > 0)
              & g["ddQTcF"].notna() & g["dQTcF"].notna()]
        c_cmax = float(np.exp(np.log(s.groupby("USUBJID")["conc"].max()).mean()))
        r_cmax = float(np.exp(np.log(o.groupby("RANDID")["conc_ngml"].max()).mean()))
        common = min(c_cmax, r_cmax)
        rd_dd = _fit(o, "ddQTcF", "conc_ngml", "RANDID", common, dd["QTcF"].provenance)
        ci_dd = _fit(s, "ddQTcF", "conc", "USUBJID", common, ddC.provenance)
        rd_d = _fit(o, "dQTcF", "conc_ngml", "RANDID", common, stepR.provenance, baseline_only=True)
        ci_d = _fit(s, "dQTcF", "conc", "USUBJID", common, stepC.provenance, baseline_only=True)
        rows.append(dict(
            drug=cipa_arm, common_conc_ngml=round(common, 3),
            ecgrdvq_ddQTcF=round(rd_dd[0], 1), ecgcipa_ddQTcF=round(ci_dd[0], 1),
            shift_ddQTcF=round(ci_dd[0] - rd_dd[0], 1),
            ecgrdvq_dQTcF_only=round(rd_d[0], 1), ecgcipa_dQTcF_only=round(ci_d[0], 1),
            shift_dQTcF_only=round(ci_d[0] - rd_d[0], 1),
            ecgrdvq_mean_placebo_applied=round(float(o["plc_applied"].mean()), 1),
            ecgcipa_mean_placebo_applied=round(float(s["plc_applied"].mean()), 1),
            diff_placebo_applied=round(float(s["plc_applied"].mean()
                                             - o["plc_applied"].mean()), 1),
            n_ecgrdvq=int(o["RANDID"].nunique()), n_ecgcipa=int(s["USUBJID"].nunique())))
    T = pd.DataFrame(rows)
    T.to_csv(out / "placebo_commonmode_test.csv", index=False)
    print(T.to_string(index=False))

    # placebo groups themselves, post-dose
    pr = g[(g["arm"] == "Placebo") & (g["TPT"] > 0)]["dQTcF"]
    pc = a[(a["TRTA"] == "Placebo") & (~a["is_baseline"])]["dQTcF"]
    pc_by_t = a[(a["TRTA"] == "Placebo") & (~a["is_baseline"])].groupby("ATPTN")["dQTcF"].agg(["mean", "count"])
    print("\nplacebo post-dose mean dQTcF  ECGRDVQ %.2f (n=%d records)   ECGCIPA %.2f (n=%d records, %d subjects)" % (
        pr.mean(), pr.size, pc.mean(), pc.size, a[a['TRTA']=='Placebo']['USUBJID'].nunique()))
    print("ECGCIPA placebo mean dQTcF by timepoint:\n", pc_by_t.round(2).T.to_string())

    # ---- J-Tpeakc sensitivity: ECGRDVQ with the power-law correction --------
    jp = jtpeakc_power(g["JTPEAK"], g["RR"])
    ddj = double_delta(g, jp, subject_col="RANDID", arm_col="arm",
                       baseline_mask=g["is_baseline"], control_arm="Placebo",
                       time_col="TPT", match_on_time=False, paired=True, name="ddJTpeakc")
    g["ddJTpeakc_pow"] = ddj.values.to_numpy()
    res = []
    for arm in ("Dofetilide", "Quinidine Sulph", "Ranolazine", "Verapamil HCL"):
        o = g[(g["arm"] == arm) & g["conc_ngml"].notna() & (g["TPT"] > 0) & g["ddJTpeakc_pow"].notna()]
        cmax, _n = cr.gmean_cmax(o)
        e = _fit(o, "ddJTpeakc_pow", "conc_ngml", "RANDID", cmax, ddj.provenance)
        res.append(dict(arm=arm, ddJTpeakc_power=round(e[0], 1), lo=round(e[1], 1), hi=round(e[2], 1)))
    J = pd.DataFrame(res); J.to_csv(out / "ecgrdvq_jtpeakc_power_sensitivity.csv", index=False)
    print("\nECGRDVQ ddJTpeakc with the power-law correction (linear-form values: 38.7, 27.7, 2.8, 0.6):")
    print(J.to_string(index=False))


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--out", type=Path, default=Path("results"))
    main(ap.parse_args().out)
