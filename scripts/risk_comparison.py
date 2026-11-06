#!/usr/bin/env python
"""TdP-risk comparison: qNet against the hiPSC-CM endpoints.

    python scripts/risk_comparison.py --wells data/cipa_myocyte_wells.csv --out results/

Reproduces results/dish_predictors_cipa12.csv, risk_classification_cipa12.csv,
risk_classification_metrics.csv, risk28_full.csv, risk28_metrics.csv and
risk28_metrics_ci.csv from the qNet outputs and the well-level table.

The 12-drug comparison is IN-SAMPLE for qNet (those are its training drugs);
the 28-drug one holds 16 out. See docs/RISK_CLASSIFICATION_28.md.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

PRED = {"qNet (in-silico)": "score_qnet",
        "arrhythmia-event fraction (dish)": "score_ead",
        "max ddFPDc (dish)": "score_maxfpd",
        "ddFPDc at Cmax (dish)": "score_cmaxfpd"}


def auc(score, pos):
    s1, s0 = score[pos], score[~pos]
    if len(s1) == 0 or len(s0) == 0:
        return np.nan
    return float((np.greater.outer(s1, s0).sum()
                  + 0.5 * np.equal.outer(s1, s0).sum()) / (len(s1) * len(s0)))


def normalise(name: pd.Series) -> pd.Series:
    """Match the wells table's drug names to the CiPA tables."""
    return (name.str.lower().str.replace(r"[^a-z]", "", regex=True)
                .replace({"dlsotalol": "sotalol"}))


def dish_predictors(wells: pd.DataFrame, drugs: pd.DataFrame) -> pd.DataFrame:
    """Blinova's predictors, re-expressed on the CiPA therapeutic-concentration axis."""
    w = wells.copy()
    w["dl"] = normalise(w["drug"])
    w = w[w.dl.isin(drugs.drug)].merge(drugs[["drug", "therapeutic"]],
                                       left_on="dl", right_on="drug", suffixes=("", "_c"))
    w["folds"] = (w["conc_um"] * 1000.0) / w["therapeutic"]
    rows = []
    for d, g in w.groupby("dl"):
        p = (g.dropna(subset=["ddFPDc"]).groupby("conc")
               .agg(m=("ddFPDc", "mean"), f=("folds", "mean")).sort_values("f"))
        rows.append(dict(
            drug=d,
            ead_fraction=g["EAD"].mean(),                       # predictors 1 and 2
            dish_max_ddFPDc=p.m.max() if len(p) else np.nan,    # predictor 4
            dish_at_cmax=float(np.interp(0.0, np.log10(p.f), p.m)) if len(p) >= 2 else np.nan,
            cmax_in_range=bool(len(p) and p.f.min() <= 1.0 <= p.f.max())))
    return pd.DataFrame(rows)


def metrics(df: pd.DataFrame, label: str) -> pd.DataFrame:
    out = []
    for lab, col in PRED.items():
        s, r = df[col].values, df["CiPA"].values
        out.append(dict(subset=label, n=len(df), predictor=lab,
                        spearman=stats.spearmanr(s, r).statistic,
                        auc_high=auc(s, r == 2), auc_low=auc(-s, r == 0)))
    return pd.DataFrame(out)


def bootstrap_ci(df, col, fn, B=4000, seed=7):
    rng = np.random.default_rng(seed)
    vals = []
    for _ in range(B):
        d = df.iloc[rng.integers(0, len(df), len(df))]
        if d["CiPA"].nunique() < 2:
            continue
        v = fn(d[col].values, d["CiPA"].values)
        if v == v:
            vals.append(v)
    return np.nanpercentile(vals, [2.5, 97.5])


def main(wells_path: Path, out: Path) -> None:
    out.mkdir(parents=True, exist_ok=True)
    wells = pd.read_csv(wells_path)

    q28 = pd.read_csv(out / "qnet28_by_drug.csv")
    D28 = dish_predictors(wells, q28[["drug", "therapeutic"]])
    X = q28.merge(D28, on="drug")
    X["score_qnet"] = -X["qnet_median"]
    X["score_ead"] = X["ead_fraction"]
    X["score_maxfpd"] = X["dish_max_ddFPDc"]
    X["score_cmaxfpd"] = X["dish_at_cmax"]
    X.to_csv(out / "risk28_full.csv", index=False)

    m = pd.concat([metrics(X, "all 28"),
                   metrics(X[X.split == "validation"], "validation 16"),
                   metrics(X[X.split == "training"], "training 12")], ignore_index=True)
    m.to_csv(out / "risk28_metrics.csv", index=False)

    ci = []
    for subset, d in [("all 28", X), ("validation 16", X[X.split == "validation"])]:
        for lab, col in PRED.items():
            e = dict(subset=subset, predictor=lab)
            for nm, fn in [("spearman", lambda s, r: stats.spearmanr(s, r).statistic),
                           ("auc_high", lambda s, r: auc(s, r == 2)),
                           ("auc_low", lambda s, r: auc(-s, r == 0))]:
                pt = fn(d[col].values, d["CiPA"].values)
                lo, hi = bootstrap_ci(d, col, fn)
                e[nm] = f"{pt:.3f} [{lo:.2f},{hi:.2f}]"
            ci.append(e)
    pd.DataFrame(ci).to_csv(out / "risk28_metrics_ci.csv", index=False)

    # the earlier 12-drug (in-sample) comparison, kept for the record
    q12 = pd.read_csv(out / "qnet_by_drug.csv")
    D12 = dish_predictors(wells, q12[["drug", "therapeutic"]])
    D12.rename(columns={"dish_at_cmax": "dish_ddFPDc_at_cmax",
                        "cmax_in_range": "cmax_in_tested_range"}).to_csv(
        out / "dish_predictors_cipa12.csv", index=False)
    X12 = q12.merge(D12, on="drug")
    X12["score_qnet"] = -X12["qnet_mean_1to4"]
    X12["score_ead"] = X12["ead_fraction"]
    X12["score_maxfpd"] = X12["dish_max_ddFPDc"]
    X12["score_cmaxfpd"] = X12["dish_at_cmax"]
    metrics(X12, "training 12").to_csv(out / "risk_classification_metrics.csv", index=False)
    print(m.to_string(index=False))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--wells", type=Path, default=Path("data/cipa_myocyte_wells.csv"))
    ap.add_argument("--out", type=Path, default=Path("results"))
    a = ap.parse_args()
    main(a.wells, a.out)
