#!/usr/bin/env python
"""Third clinical corpus: the CiPA ECG validation study, and cross-study replication.

    python scripts/ecgcipa_reference.py --out results/

Writes results/ecgcipa_concentration_qtc.csv,
results/ecgcipa_correction_validation.csv and
results/cross_study_replication.csv.

Why this study is here. ECGRDVQ and ECGDMMLD between them cover 8 molecules in 44
healthy volunteers, and every drug in the dish-to-clinic bridge rests on a single
study. PhysioNet ECGCIPA (Vicente J et al., Clin Pharmacol Ther 2019;105:943-953,
doi:10.1002/cpt.1303) adds 60 subjects and three new molecules, and - more useful -
gives ranolazine, verapamil and dofetilide a second, independent estimate from a
different design.

It also ships the study's **own** derived QTcF and J-Tpeakc per record, plus
change-from-baseline, which makes it a per-record check on the correction code.

Design note: five arms (ranolazine, verapamil, lopinavir+ritonavir, chloroquine,
placebo) are parallel with 10 subjects each; dofetilide and diltiazem+dofetilide
form a 2-period crossover in a further 10, with no placebo period of its own.
The single part 1 placebo group is subtracted from every arm, so placebo
subtraction is **unpaired for all of them** and its sampling error is shared by
all ECGCIPA drug effects. See scripts/placebo_commonmode_test.py for the measured
size of that effect (about 2-4 ms of a 5-10 ms cross-study shift).
"""
from __future__ import annotations

import argparse
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

from dish2ecg import (fit_concentration_response, fridericia, jtpeakc_power,
                      subtract_baseline, subtract_control)
from dish2ecg.io import ecgcipa_wide, load_ecgcipa
from dish2ecg.provenance import Endpoint

#: arm -> the analyte whose concentration drives that arm's exposure-response
#: Design of each ECGCIPA arm, and how placebo is subtracted. There is ONE placebo
#: group (part 1, 10 subjects) and it is subtracted, time-matched, from every
#: arm - including the part 2 dofetilide crossover arm, which has no placebo
#: period of its own. So placebo subtraction is unpaired for all of them.
ECGCIPA_DESIGN = {
    "Ranolazine": "parallel, 10/arm (unpaired placebo)",
    "Verapamil": "parallel, 10/arm (unpaired placebo)",
    "Dofetilide": "part 2 crossover, 10 subjects (placebo = part 1 group, unpaired)",
}


ARM_ANALYTE = {
    "Ranolazine": "Ranolazine", "Verapamil": "Verapamil",
    "Chloroquine": "Chloroquine", "Lopinavir+Ritonavir": "Lopinavir",
    "Dofetilide": "Dofetilide", "Diltiazem+Dofetilide": "Dofetilide",
}
#: ECGCIPA arm -> ECGRDVQ arm, for the three drugs present in both studies
SHARED = {"Ranolazine": "Ranolazine", "Verapamil": "Verapamil HCL",
          "Dofetilide": "Dofetilide"}
_UNIT = {"pg/ml": 1e-3, "ng/ml": 1.0, "ug/ml": 1e3, "\u00b5g/ml": 1e3}


def analysis_set(adeg: pd.DataFrame) -> pd.DataFrame:
    """One row per subject x period x timepoint, using the study's triplicate mean.

    ADEG carries the three replicate ECGs (``EGREPNUM`` 1-3) *and* the study's own
    per-timepoint average of them (``EGREPNUM`` 0 after filling). Taking row 0
    rather than re-averaging keeps the analysis identical to the study's own.
    """
    w = ecgcipa_wide(adeg).dropna(subset=["QT", "RR", "JTP"]).copy()
    a = w[w["EGREPNUM"] == 0].reset_index(drop=True)
    a["is_baseline"] = a["ABLFL"].eq("Y")
    return a


def concentrations_ngml(adpc: pd.DataFrame, analyte: str) -> pd.DataFrame:
    """Plasma concentration of one analyte, harmonised to ng/mL.

    Units differ **by analyte** in this study: dofetilide is pg/mL, every other
    analyte ng/mL. Comparing across studies without this conversion silently
    mixes scales by a factor of 1000.
    """
    m = adpc.rename(columns={"AVAL": "conc", "PARAM": "analyte",
                             "PCSTRESU": "unit"})
    m = m[m["analyte"] == analyte].copy()
    if m.empty:
        raise KeyError(f"analyte {analyte!r} not in ADPC")
    u = str(m["unit"].iloc[0]).lower()
    if u not in _UNIT:
        raise ValueError(f"unhandled concentration unit {u!r} for {analyte}")
    m["conc"] = m["conc"] * _UNIT[u]
    out = m.groupby(["USUBJID", "ATPTN"], as_index=False)["conc"].mean()
    out.attrs["source_unit"] = u
    return out


def validate_against_study_values(a: pd.DataFrame, adeg: pd.DataFrame) -> pd.DataFrame:
    """Check our QTcF, J-Tpeakc and change-from-baseline against the study's own."""
    rows = []
    ours_qtcf = fridericia(a["QT"], a["RR"]).values.to_numpy()
    ours_jtpc = jtpeakc_power(a["JTP"], a["RR"]).values.to_numpy()
    for label, ours, theirs in [("QTcF", ours_qtcf, a["QTCF"].to_numpy()),
                                ("JTpeakc", ours_jtpc, a["JTPC"].to_numpy())]:
        ok = np.isfinite(ours) & np.isfinite(theirs)
        e = ours[ok] - theirs[ok]
        rows.append(dict(quantity=label, n=int(ok.sum()),
                         median_diff_ms=round(float(np.median(e)), 4),
                         sd_ms=round(float(np.std(e, ddof=1)), 4),
                         max_abs_ms=round(float(np.max(np.abs(e))), 3)))
    # the exponent their J-Tpeakc implies
    ok = (a["JTP"] > 0) & (a["JTPC"] > 0) & (a["RR"] > 0)
    alpha = float(np.polyfit(np.log(a.loc[ok, "RR"] / 1000.0),
                             np.log(a.loc[ok, "JTP"] / a.loc[ok, "JTPC"]), 1)[0])
    rows.append(dict(quantity="JTpeak exponent implied by study", n=int(ok.sum()),
                     median_diff_ms=round(alpha, 4), sd_ms=np.nan, max_abs_ms=np.nan))
    return pd.DataFrame(rows)


def main(args) -> None:
    warnings.filterwarnings("ignore")
    args.out.mkdir(parents=True, exist_ok=True)
    D = load_ecgcipa()
    a = analysis_set(D["adeg"])
    print(f"ECGCIPA: {len(a)} analysis records, {a.USUBJID.nunique()} subjects, "
          f"{a.TRTA.nunique()} arms")

    V = validate_against_study_values(a, D["adeg"])
    V.to_csv(args.out / "ecgcipa_correction_validation.csv", index=False)
    print("\n" + V.to_string(index=False))

    qtcf = fridericia(a["QT"], a["RR"])
    jtpc = jtpeakc_power(a["JTP"], a["RR"])
    dd = {}
    for nm, ep in (("QTcF", qtcf), ("JTpeakc", jtpc)):
        step = subtract_baseline(a, ep, subject_col="USUBJID", arm_col="TRTA",
                                 baseline_mask=a["is_baseline"], match_on_time=False)
        dd[nm] = subtract_control(a, step, arm_col="TRTA", control_arm="Placebo",
                                  time_col="ATPTN", paired=False, name="dd" + nm)
        a["dd" + nm] = dd[nm].values.to_numpy()

    rows = []
    for arm, analyte in ARM_ANALYTE.items():
        pk = concentrations_ngml(D["adpc"], analyte)
        s = a[(a["TRTA"] == arm) & (~a["is_baseline"])].merge(
            pk, on=["USUBJID", "ATPTN"], how="inner")
        s = s[s["conc"].notna()]
        if len(s) < 10:
            continue
        cmax = float(np.exp(np.log(s.groupby("USUBJID")["conc"].max()).mean()))
        for nm in ("QTcF", "JTpeakc"):
            sub = s[s["dd" + nm].notna()]
            ep = Endpoint(values=sub["dd" + nm].astype(float), name="dd" + nm,
                          units="ms", domain="clinical", provenance=dd[nm].provenance)
            f = fit_concentration_response(sub, ep, concentration_col="conc",
                                           group_col="USUBJID", at_concentration=cmax)
            r = f.to_row()
            r.update(study="ecgcipa", arm=arm, analyte=analyte,
                     gmean_cmax_ngml=cmax, conc_units="ng/mL")
            rows.append(r)
    E = pd.DataFrame(rows)
    E["exceeds_10ms"] = E["effect_ci90_upper"] > 10.0
    E.to_csv(args.out / "ecgcipa_concentration_qtc.csv", index=False)
    print("\n" + E[["arm", "endpoint", "n_obs", "n_groups", "gmean_cmax_ngml",
                     "effect_at_concentration", "effect_ci90_lower",
                     "effect_ci90_upper"]].round(2).to_string(index=False))

    # cross-study replication at a matched plasma concentration
    clin = pd.read_csv(args.clinical)
    rd = clin[clin["study"] == "ecgrdvq"]
    rows = []
    for cipa_arm, rdvq_arm in SHARED.items():
        pk = concentrations_ngml(D["adpc"], ARM_ANALYTE[cipa_arm])
        s = a[(a["TRTA"] == cipa_arm) & (~a["is_baseline"])].merge(
            pk, on=["USUBJID", "ATPTN"], how="inner")
        s = s[s["conc"].notna() & s["ddQTcF"].notna()]
        o = rd[(rd["arm"] == rdvq_arm) & rd["conc_ngml"].notna() & (rd["TPT"] > 0)
               & rd["ddQTcF"].notna()]
        c_cmax = float(np.exp(np.log(s.groupby("USUBJID")["conc"].max()).mean()))
        r_cmax = float(np.exp(np.log(o.groupby("RANDID")["conc_ngml"].max()).mean()))
        common = min(c_cmax, r_cmax)

        def _fit(frame, ccol, gcol):
            ep = Endpoint(values=frame["ddQTcF"].astype(float), name="ddQTcF",
                          units="ms", domain="clinical",
                          provenance=dd["QTcF"].provenance)
            return fit_concentration_response(frame, ep, concentration_col=ccol,
                                              group_col=gcol,
                                              at_concentration=common)

        fc, fr = _fit(s, "conc", "USUBJID"), _fit(o, "conc_ngml", "RANDID")
        rows.append(dict(
            drug=cipa_arm, common_conc_ngml=round(common, 3),
            ecgrdvq_design="crossover (paired placebo)",
            ecgcipa_design=ECGCIPA_DESIGN[cipa_arm],
            ecgrdvq_n=fr.n_groups, ecgcipa_n=fc.n_groups,
            ecgrdvq_ddQTcF=round(fr.effect_at_concentration, 1),
            ecgrdvq_lo=round(fr.effect_ci90[0], 1), ecgrdvq_hi=round(fr.effect_ci90[1], 1),
            ecgcipa_ddQTcF=round(fc.effect_at_concentration, 1),
            ecgcipa_lo=round(fc.effect_ci90[0], 1), ecgcipa_hi=round(fc.effect_ci90[1], 1),
            difference_ms=round(fc.effect_at_concentration
                                - fr.effect_at_concentration, 1),
            intervals_overlap=bool(fc.effect_ci90[0] <= fr.effect_ci90[1]
                                   and fr.effect_ci90[0] <= fc.effect_ci90[1])))
    X = pd.DataFrame(rows)
    X.to_csv(args.out / "cross_study_replication.csv", index=False)
    print("\n" + X.to_string(index=False))
    print(f"\nbetween-study difference at matched concentration: "
          f"{X.difference_ms.min():+.1f} to {X.difference_ms.max():+.1f} ms; "
          f"intervals overlap for {int(X.intervals_overlap.sum())} of {len(X)} drugs")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=Path("results"))
    ap.add_argument("--clinical", type=Path,
                    default=Path("results/clinical_ddqtc.csv"))
    main(ap.parse_args())
