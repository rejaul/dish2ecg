#!/usr/bin/env python
"""Does the T-end rule explain the ECGRDVQ-vs-ECGCIPA shift?  (It does not.)

    python scripts/delineation_shift.py --wf wf2/ --out results/

`--wf` holds median-beat records as <study>/<EGREFID>.{hea,dat,atr} for
study in {ecgrdvq, ecgcipa}; fetch them from PhysioNet (nothing is redistributed).

Writes results/delineation_both.csv and results/delineation_cascade_effect.csv.

The baseline/placebo cascade is LINEAR in QT, so the amount a published ddQTcF
would move if the T-end rule changed is exactly the cascade applied to the
per-record difference QTcF_tangent - QTcF_study. That is what this computes.
"""
from __future__ import annotations
import argparse, warnings
from pathlib import Path
import numpy as np
import pandas as pd

from dish2ecg.viz import read_wfdb, read_wfdb_annotations
from dish2ecg.io import load_physionet_clinical, load_ecgcipa

ARMS = ("Ranolazine", "Dofetilide")
PLACEBO = "Placebo"


def fridericia_qtc(qt, rr):
    return qt / (rr / 1000.0) ** (1 / 3)


def tangent_tend(y, t_peak_ms, fs):
    """Tangent (Lepeschkin) T-end: steepest T downslope extrapolated to baseline."""
    i0 = int(t_peak_ms * fs / 1000.0)
    i1 = min(len(y) - 2, i0 + int(0.42 * fs))
    if i1 - i0 < 10:
        return np.nan
    seg = y[i0:i1]
    d = np.gradient(seg)
    k = int(np.argmin(d))
    if d[k] >= -1e-9:
        return np.nan
    slope = d[k] * fs / 1000.0
    base = np.median(y[max(0, i0 - int(0.38 * fs)):max(1, i0 - int(0.30 * fs))])
    return (i0 + k) * 1000.0 / fs + (base - seg[k]) / slope


def fiducials(ann, fs):
    ms = [(t * 1000.0 / fs, c) for t, c in ann]
    on = [t for t, c in ms if c == "("]
    off = [t for t, c in ms if c == ")"]
    tpk = [t for t, c in ms if c == "t"]
    if len(on) < 2 or len(off) < 2 or not tpk:
        return None
    return dict(QRS_onset=on[1], T_peak=tpk[0], T_end=off[1])


def delineate(wf: Path) -> pd.DataFrame:
    rows = []
    for study in ("ecgrdvq", "ecgcipa"):
        for f in sorted((wf / study).glob("*.atr")):
            try:
                rec = read_wfdb(f.with_suffix(""))
                F = fiducials(read_wfdb_annotations(f), rec.fs)
                if F is None:
                    continue
                te = tangent_tend(rec.lead("VCGMAG"), F["T_peak"], rec.fs)
                rows.append(dict(study=study, EGREFID=f.stem, qrs_onset=F["QRS_onset"],
                                 t_peak=F["T_peak"], qt_study=F["T_end"] - F["QRS_onset"],
                                 qt_tangent=(te - F["QRS_onset"]) if te == te else np.nan))
            except Exception:
                continue
    D = pd.DataFrame(rows)
    D["dQT"] = D.qt_tangent - D.qt_study
    D["plausible"] = D.qt_tangent.between(250, 650)
    return D


def double_delta_delta(d: pd.DataFrame, drug: str, B: int = 2000, seed: int = 11):
    """Baseline- then placebo-corrected delineation difference, per subject."""
    per = {}
    for arm in (drug, PLACEBO):
        a = d[d.arm == arm]
        if a.empty:
            return None
        bl = a[a.base].groupby("subj").delta.median()
        po = a[~a.base].groupby("subj").delta.median()
        common = po.index.intersection(bl.index)
        if not len(common):
            return None
        per[arm] = (po.loc[common] - bl.loc[common]).values
    rng = np.random.default_rng(seed)
    dr, pl = per[drug], per[PLACEBO]
    bs = [np.median(rng.choice(dr, len(dr))) - np.median(rng.choice(pl, len(pl))) for _ in range(B)]
    lo, hi = np.percentile(bs, [2.5, 97.5])
    return dict(drug=drug, n_drug=len(dr), n_plac=len(pl),
                ddelta=float(np.median(dr) - np.median(pl)), lo=float(lo), hi=float(hi))


def main(wf: Path, out: Path) -> None:
    warnings.filterwarnings("ignore")
    out.mkdir(parents=True, exist_ok=True)
    D = delineate(wf)
    D.to_csv(out / "delineation_both.csv", index=False)

    rd = load_physionet_clinical("ecgrdvq")
    r1 = rd[["EGREFID", "RANDID", "EXTRT", "BASELINE", "RR"]].merge(
        D[D.study == "ecgrdvq"][["EGREFID", "qt_study", "qt_tangent"]], on="EGREFID")
    r1 = r1[r1.qt_tangent.between(250, 650)].rename(columns={"RANDID": "subj", "EXTRT": "arm"})
    r1["base"] = r1.BASELINE.astype(str).str.upper().eq("Y")
    r1["delta"] = fridericia_qtc(r1.qt_tangent, r1.RR) - fridericia_qtc(r1.qt_study, r1.RR)

    eg = load_ecgcipa()["adeg"]
    need = eg[eg.PARAMCD.isin(["QT", "RR"])][["EGREFID", "USUBJID", "TRTA", "ATPTN", "PARAMCD", "AVAL"]]
    w2 = need.pivot_table(index=["EGREFID", "USUBJID", "TRTA", "ATPTN"], columns="PARAMCD",
                          values="AVAL", aggfunc="mean").reset_index().merge(
        D[D.study == "ecgcipa"][["EGREFID", "qt_study", "qt_tangent"]], on="EGREFID")
    w2 = w2[w2.qt_tangent.between(250, 650) & w2.RR.notna()].rename(
        columns={"USUBJID": "subj", "TRTA": "arm"})
    # ECGCIPA's ABLFL='Y' row is a DERIVED pre-dose average with no waveform;
    # the real baseline records are the three pre-dose timepoints (-1, -0.5, 0 h).
    w2["base"] = w2.ATPTN <= 3
    w2["delta"] = fridericia_qtc(w2.qt_tangent, w2.RR) - fridericia_qtc(w2.qt_study, w2.RR)

    rows = [dict(study=nm, **r)
            for nm, d in (("ecgrdvq", r1), ("ecgcipa", w2))
            for a in ARMS
            if (r := double_delta_delta(d, a)) is not None]
    DD = pd.DataFrame(rows)
    DD.to_csv(out / "delineation_cascade_effect.csv", index=False)
    print(DD.to_string(index=False))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--wf", type=Path, required=True)
    ap.add_argument("--out", type=Path, default=Path("results"))
    a = ap.parse_args()
    main(a.wf, a.out)
