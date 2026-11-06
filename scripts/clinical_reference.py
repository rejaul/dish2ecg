#!/usr/bin/env python
"""Clinical reference analysis: ECGRDVQ + ECGDMMLD -> fully corrected endpoints.

Reproduces results/clinical_ddqtc.csv, results/concentration_qtc_results.csv,
results/rate_correction_exponents.csv and results/validation_vs_published.csv.

    python scripts/clinical_reference.py --out results/

Data are downloaded from PhysioNet on first run and cached; nothing is
redistributed with this repository.
"""
from __future__ import annotations

import argparse
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from dish2ecg import (bazett, double_delta, fit_concentration_response,
                      fit_individual_exponents, fit_power_exponent, fridericia,
                      individual_correction, jtpeakc_hnatkova2019)
from dish2ecg.io import load_physionet_clinical
from dish2ecg.provenance import CorrectionRecord, Endpoint

# Baseline- and placebo-corrected change at the population Cmax reported for
# ECGRDVQ in Vicente J et al., J Am Heart Assoc 2015;4:e001615, Table 2
# (doi:10.1161/JAHA.114.001615). None = the paper reports no concentration
# relationship for that biomarker.
PUBLISHED_RDVQ = {
    ("Dofetilide", "ddQTcF"): (73.6, 65.8, 81.5),
    ("Quinidine Sulph", "ddQTcF"): (78.9, 68.2, 89.7),
    ("Ranolazine", "ddQTcF"): (12.0, 7.3, 16.7),
    ("Verapamil HCL", "ddQTcF"): None,
    ("Dofetilide", "ddJTpeakc"): (39.1, 31.6, 46.6),
    ("Quinidine Sulph", "ddJTpeakc"): (26.1, 13.5, 38.7),
    ("Ranolazine", "ddJTpeakc"): None,
    ("Verapamil HCL", "ddJTpeakc"): None,
}

# Placebo- and baseline-corrected changes at Tmax reported for ECGDMMLD in
# Vicente J et al., PLoS ONE 2016;11:e0163619, Table 3 (mean [95% CI], ms).
PUBLISHED = {
    ("Dofetilide", "ddQTcF"): (45.2, 37.8, 52.6),
    ("Dofetilide", "ddJTpeakc"): (29.1, 20.2, 38.0),
    ("Dofetilide", "ddTpeakTend"): (17.9, 9.0, 26.7),
    ("Mexiletine + Dofetilide", "ddQTcF"): (17.3, 11.6, 23.0),
    ("Mexiletine + Dofetilide", "ddJTpeakc"): (1.1, -5.5, 7.7),
    ("Mexiletine + Dofetilide", "ddTpeakTend"): (14.6, 10.0, 19.2),
}

# Wide concentration columns of ECGDMMLD: arm -> (column, analyte, unit)
DMMLD_CONC = {
    "Dofetilide": ("DOF", "Dofetilide", "pg/mL"),
    "Lidocaine + Dofetilide": ("DOF", "Dofetilide", "pg/mL"),
    "Mexiletine + Dofetilide": ("DOF", "Dofetilide", "pg/mL"),
    "Moxifloxacin + Diltiazem": ("MOXI", "Moxifloxacin", "ng/mL"),
}

KEY = ["study", "RANDID", "VISIT", "arm", "TPT", "is_baseline"]


def prepare(df: pd.DataFrame, study: str) -> tuple[pd.DataFrame, dict]:
    """Repair known defects and attach a single concentration column per row."""
    d = df.copy()
    d["study"] = study
    d["arm"] = d[df.attrs["arm_col"]]

    # Two ECGRDVQ PR values are stored as 32-bit wrapped negatives. An unguarded
    # mean of PR is about -1.6 million ms; adding 2**32 recovers 345 and 324 ms.
    n_wrap = int((d["PR"] < -1000).sum())
    d.loc[d["PR"] < -1000, "PR"] += 2 ** 32

    # QT == QRS + JTPEAK + TPEAKTEND holds for every complete row in both studies,
    # so JTPEAK is reconstructable where it is missing.
    m = d["JTPEAK"].isna() & d[["QT", "QRS", "TPEAKTEND"]].notna().all(axis=1)
    n_rec = int(m.sum())
    d.loc[m, "JTPEAK"] = d.loc[m, "QT"] - d.loc[m, "QRS"] - d.loc[m, "TPEAKTEND"]

    d["is_baseline"] = d["BASELINE"].astype(str).str.upper().eq("Y")

    if df.attrs["concentration_layout"] == "long":
        d["conc_raw"] = d["PCSTRESN"]
        d["analyte"] = d["PCTEST"].astype("object")
        d["conc_unit"] = d["PCSTRESU"].astype("object")
    else:
        d["conc_raw"] = np.nan
        d["analyte"] = pd.Series([None] * len(d), dtype="object")
        d["conc_unit"] = pd.Series([None] * len(d), dtype="object")
        for arm, (col, analyte, unit) in DMMLD_CONC.items():
            sel = d["arm"].eq(arm).to_numpy()
            d.loc[sel, "conc_raw"] = pd.to_numeric(d.loc[sel, col], errors="coerce")
            d.loc[sel, "analyte"] = analyte
            d.loc[sel, "conc_unit"] = unit

    # Units are not uniform: dofetilide is pg/mL in both studies, everything else
    # ng/mL. Harmonise before any pooled calculation.
    u = d["conc_unit"].astype(str).str.lower()
    conc = d["conc_raw"].astype(float).copy()
    conc[u.eq("pg/ml")] = conc[u.eq("pg/ml")] / 1000.0
    d["conc_ngml"] = conc

    return d, {"pr_wraps_repaired": n_wrap, "jtpeak_reconstructed": n_rec}


def collapse_replicates(d: pd.DataFrame) -> pd.DataFrame:
    """Average the three replicate ECGs per subject x visit x timepoint (ICH E14)."""
    num = ["RR", "QT", "QRS", "JTPEAK", "TPEAKTEND", "conc_ngml"]
    g = d.groupby(KEY, dropna=False, as_index=False).agg(
        **{c: (c, "mean") for c in num}, n_rep=("EGREFID", "size"))
    g["analyte"] = d.groupby(KEY, dropna=False)["analyte"].first().to_numpy()
    return g


def corrected_endpoints(g: pd.DataFrame) -> tuple[dict, dict, pd.DataFrame, float]:
    """All rate-corrected endpoints, then the baseline/placebo cascade."""
    drugfree = g["is_baseline"] | g["arm"].eq("Placebo")
    alpha_pop, diag = fit_power_exponent(g.loc[drugfree, "QT"], g.loc[drugfree, "RR"])
    ex = fit_individual_exponents(g[drugfree], interval_col="QT",
                                  cycle_length_col="RR", subject_col="RANDID",
                                  min_records=10)
    eps = {
        "QTcF": fridericia(g["QT"], g["RR"]),
        "QTcB": bazett(g["QT"], g["RR"], name="QTcB"),
        "QTcI": individual_correction(g, interval_col="QT", cycle_length_col="RR",
                                      subject_col="RANDID", exponents=ex,
                                      fallback_alpha=alpha_pop, name="QTcI"),
        "JTpeakc": jtpeakc_hnatkova2019(g["JTPEAK"], g["RR"]),
        # Tpeak-Tend was not rate-corrected in the published analysis.
        "TpeakTend": Endpoint(
            values=g["TPEAKTEND"].astype(float).rename("TpeakTend"),
            name="TpeakTend", units="ms", domain="clinical",
            provenance=(CorrectionRecord(
                step="rate_correction", method="none",
                formula="no rate correction applied",
                note="Tpeak-Tend was not rate-corrected in Vicente 2016."),)),
    }
    dd = {nm: double_delta(g, ep, subject_col="RANDID", arm_col="arm",
                           baseline_mask=g["is_baseline"], control_arm="Placebo",
                           time_col="TPT",
                           # one pre-dose baseline per subject-period, not a
                           # day-matched baseline across periods
                           match_on_time=False, paired=True, name="dd" + nm)
          for nm, ep in eps.items()}
    return eps, dd, ex, alpha_pop


def gmean_cmax(d: pd.DataFrame) -> tuple[float, int]:
    per = d.groupby("RANDID")["conc_ngml"].max().dropna()
    per = per[per > 0]
    return float(np.exp(np.log(per).mean())), int(per.size)


def main(out: Path) -> None:
    warnings.filterwarnings("ignore")
    out.mkdir(parents=True, exist_ok=True)

    frames, defects = [], {}
    for study in ("ecgrdvq", "ecgdmmld"):
        raw = load_physionet_clinical(study)
        d, defect = prepare(raw, study)
        defects[study] = defect
        frames.append(collapse_replicates(d))
    A = pd.concat(frames, ignore_index=True)
    print("data defects repaired:", defects)

    per_study, fits, exps = {}, [], []
    for study, g0 in A.groupby("study"):
        g = g0.reset_index(drop=True)
        eps, dd, ex, alpha_pop = corrected_endpoints(g)
        per_study[study] = (g, eps, dd)

        e = ex.copy(); e["study"] = study; e["pooled_alpha"] = alpha_pop
        exps.append(e)

        for arm in sorted(g["arm"].unique()):
            if arm == "Placebo":
                continue
            sel = ((g["arm"] == arm) & g["conc_ngml"].notna()
                   & (g["TPT"] > 0)).to_numpy()
            if sel.sum() < 10:
                continue
            cmax, nsub = gmean_cmax(g[sel])
            for nm in ("QTcF", "QTcI", "JTpeakc"):
                f = fit_concentration_response(
                    g[sel], dd[nm].subset(sel), concentration_col="conc_ngml",
                    group_col="RANDID", at_concentration=cmax)
                row = f.to_row()
                row.update(study=study, arm=arm, conc_units="ng/mL",
                           analyte=g.loc[sel, "analyte"].iloc[0],
                           gmean_cmax_ngml=cmax, n_subjects_cmax=nsub)
                fits.append(row)

    # long endpoint table
    longs = []
    for study, (g, eps, dd) in per_study.items():
        t = g.copy()
        for nm, ep in eps.items():
            t[nm] = ep.values.to_numpy()
        for nm, ep in dd.items():
            t["dd" + nm] = ep.values.to_numpy()
        longs.append(t)
    pd.concat(longs, ignore_index=True).to_csv(out / "clinical_ddqtc.csv", index=False)

    F = pd.DataFrame(fits)
    F["exceeds_10ms"] = F["effect_ci90_upper"] > 10.0
    F.to_csv(out / "concentration_qtc_results.csv", index=False)
    pd.concat(exps, ignore_index=True).to_csv(
        out / "rate_correction_exponents.csv", index=False)

    # validation against the published ECGDMMLD values
    g, _, dd = per_study["ecgdmmld"]
    t = g.copy()
    for nm, ep in dd.items():
        t["dd" + nm] = ep.values.to_numpy()
    tmax = (t[t["conc_ngml"].notna()].groupby(["arm", "TPT"])["conc_ngml"].mean()
            .reset_index().sort_values("conc_ngml").groupby("arm").tail(1)
            .set_index("arm")["TPT"].to_dict())
    rows = []
    for (arm, col), (pm, plo, phi) in PUBLISHED.items():
        x = t.loc[(t["arm"] == arm) & (t["TPT"] == tmax[arm]), "dd" + col[2:]]
        x = x.dropna().to_numpy()
        m = x.mean(); se = x.std(ddof=1) / np.sqrt(x.size)
        h = stats.t.ppf(0.975, x.size - 1) * se
        rows.append(dict(arm=arm, endpoint=col, tmax_h=tmax[arm], n_subjects=x.size,
                         computed=round(m, 1), computed_lo=round(m - h, 1),
                         computed_hi=round(m + h, 1), published=pm,
                         published_lo=plo, published_hi=phi, diff=round(m - pm, 1),
                         published_in_our_CI=bool(m - h <= pm <= m + h)))
    V = pd.DataFrame(rows)
    V.to_csv(out / "validation_vs_published.csv", index=False)

    # ECGRDVQ: our concentration-response fits against the published Table 2
    idx = F.set_index(["arm", "endpoint"])
    rdvq = []
    for (arm, ep), pub in PUBLISHED_RDVQ.items():
        r = idx.loc[(arm, ep)]
        if isinstance(r, pd.DataFrame):
            r = r[r["study"] == "ecgrdvq"].iloc[0]
        lo, hi = r["effect_ci90_lower"], r["effect_ci90_upper"]
        if pub is None:
            agree = bool(lo <= 0.0 <= hi)
            rdvq.append(dict(arm=arm, endpoint=ep, published="no relationship",
                             ours=round(r["effect_at_concentration"], 1),
                             ours_lo=round(lo, 1), ours_hi=round(hi, 1),
                             agreement="our 90% CI includes 0" if agree
                             else "DISAGREES", agrees=agree))
        else:
            pm, plo, phi = pub
            agree = bool(lo <= pm <= hi)
            rdvq.append(dict(arm=arm, endpoint=ep, published=pm,
                             published_lo=plo, published_hi=phi,
                             ours=round(r["effect_at_concentration"], 1),
                             ours_lo=round(lo, 1), ours_hi=round(hi, 1),
                             agreement="published point inside our 90% CI" if agree
                             else "OUTSIDE our CI", agrees=agree))
    RV = pd.DataFrame(rdvq)
    RV.to_csv(out / "validation_ecgrdvq_vs_published.csv", index=False)
    print("\n--- ECGRDVQ vs Vicente 2015 JAHA Table 2 ---")
    print(RV.to_string(index=False))
    print(f"agreement: {int(RV.agrees.sum())} of {len(RV)}")
    print(V.to_string(index=False))
    print(f"\npublished value inside our 95% CI for "
          f"{int(V.published_in_our_CI.sum())} of {len(V)} comparisons")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=Path("results"))
    main(ap.parse_args().out)
