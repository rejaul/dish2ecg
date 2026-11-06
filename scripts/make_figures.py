"""Build the explanatory figures. Real data where it exists, schematic where it does not."""
from __future__ import annotations

import argparse, warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from dish2ecg import fridericia, fit_power_exponent
from dish2ecg.io import load_physionet_clinical
from dish2ecg.viz import (read_wfdb, read_wfdb_annotations, plot_ecg_beat,
                          plot_field_potential, plot_qt_rr,
                          plot_concentration_response)

#: One median beat whose annotation intervals reproduce the clinical table exactly.
DEMO_RECORD = "00ed2097-cd14-4f03-ab33-853da5be5550"


def fiducials_from_annotations(anns, fs):
    """P onset, QRS onset, J point, T peak, T end from WFDB waveform annotations."""
    ms = [(t * 1000.0 / fs, c) for t, c in anns]
    on = [t for t, c in ms if c == "("]
    off = [t for t, c in ms if c == ")"]
    tpk = [t for t, c in ms if c == "t"]
    out = {}
    if len(on) >= 2:
        out["P_onset"], out["QRS_onset"] = on[0], on[1]
    if off:
        out["J"] = off[0]
    if tpk:
        out["T_peak"] = tpk[0]
    if len(off) >= 2:
        out["T_end"] = off[1]
    return out


def main(out: Path, wf: Path) -> None:
    warnings.filterwarnings("ignore")
    out.mkdir(parents=True, exist_ok=True)
    clin = load_physionet_clinical("ecgrdvq")

    # ---- 1. real annotated beat + schematic field potential -----------------
    fig, axes = plt.subplots(1, 2, figsize=(10.4, 3.5))
    rec = read_wfdb(wf / "beat")
    anns = read_wfdb_annotations(wf / "beat.atr")
    fid = fiducials_from_annotations(anns, rec.fs)
    row = clin[clin["EGREFID"] == DEMO_RECORD].iloc[0]
    iv = {k: float(row[k]) for k in ("PR", "QRS", "QT", "JTPEAK", "TPEAKTEND")
          if k in row and pd.notna(row[k])}
    plot_ecg_beat(rec, lead="II", fiducials=fid, intervals=iv, ax=axes[0],
                  window_ms=(0.0, 820.0), show_intervals=("PR","QRS","JTPEAK","TPEAKTEND","QT"),
                  title=(f"Real median beat - subject {int(row['RANDID'])}, "
                         f"{row['EXTRT']}, {rec.fs:g} Hz"))
    plot_field_potential(ax=axes[1])
    fig.tight_layout(); fig.savefig(out / "viz_beat_and_fp.png", dpi=200)
    plt.close(fig)
    check = iv.get("QRS", 0) + iv.get("JTPEAK", 0) + iv.get("TPEAKTEND", 0)
    print(f"beat: fiducials {({k: round(v) for k, v in fid.items()})}")
    print(f"      intervals {iv}  | QRS+JTp+TpTe = {check:.0f} vs QT {iv.get('QT'):.0f}")
    print(f"      ADC: {rec.adc_resolution}-bit, lead II gain {rec.gain[1]:.1f} "
          f"units/mV -> {rec.quantisation_uv(1)*1000:.3f} nV per LSB")

    # ---- 2. QT against RR, drug-free records, with correction curves --------
    g = clin.copy()
    base = g[(g["BASELINE"] == "Y") | (g["EXTRT"] == "Placebo")]
    base = base[np.isfinite(base["QT"]) & np.isfinite(base["RR"])]
    alpha, _ = fit_power_exponent(base["QT"], base["RR"])
    fig, ax = plt.subplots(figsize=(5.6, 3.6))
    plot_qt_rr(base["RR"], base["QT"], fitted_alpha=alpha, ax=ax,
               title="QT against RR, drug-free records (ECGRDVQ)")
    fig.tight_layout(); fig.savefig(out / "viz_qt_rr.png", dpi=200); plt.close(fig)
    print(f"QT/RR: n={len(base)}, fitted alpha={alpha:.4f}")

    # ---- 3. concentration-response for the four bridge drugs ----------------
    res = pd.read_csv(out / "concentration_qtc_results.csv")
    res = res[(res.study == "ecgrdvq") & (res.endpoint == "ddQTcF")]
    arms = ["Dofetilide", "Quinidine Sulph", "Ranolazine", "Verapamil HCL"]
    dd = pd.read_csv(out / "clinical_ddqtc.csv")
    dd = dd[dd["study"] == "ecgrdvq"]
    fig, axes = plt.subplots(1, 4, figsize=(13.2, 3.2), sharey=True)
    for ax, arm in zip(axes, arms):
        r = res[res.arm == arm]
        s = dd[(dd["arm"] == arm) & dd["conc_ngml"].notna() & dd["ddQTcF"].notna()]
        if r.empty or s.empty:
            ax.set_visible(False); continue
        r = r.iloc[0]
        plot_concentration_response(
            s["conc_ngml"], s["ddQTcF"], slope=r.slope, intercept=r.effect_at_concentration,
            at_concentration=r.gmean_cmax_ngml, effect_at=r.effect_at_concentration,
            ci=(r.effect_ci90_lower, r.effect_ci90_upper), ax=ax,
            xlabel="plasma conc (ng/mL)", title=arm,
            ylabel="$\\Delta\\Delta$QTcF (ms)" if arm == arms[0] else "")
    fig.tight_layout(); fig.savefig(out / "viz_concentration_response.png", dpi=200)
    plt.close(fig)
    print("figures written to", out)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=Path("results"))
    ap.add_argument("--wf", type=Path, default=Path("../wf"))
    main(ap.parse_args().out, ap.parse_args().wf)
