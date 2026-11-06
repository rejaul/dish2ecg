"""Two questions that only the PhysioNet WAVEFORMS can answer.

ECGRDVQ and ECGDMMLD ship waveforms as well as interval tables: `medians/` holds
one median beat per record (16 leads, 1000 Hz) with a `.atr` carrying the study's
own fiducial points, and `raw/` holds the 12-lead, 10-second strip the median was
built from. The clinical CSV column ``EGREFID`` is the join key.

**1. Can QT/RR hysteresis be corrected here?** No, and the reason is structural
rather than a missing implementation. QT needs one to two minutes to settle after
a rate change, so the RR history that sets the QT in a 10-second strip mostly
precedes the strip. This script quantifies what that costs: how much RR moves
inside the window, and the QTcF spread that the choice of pairing RR implies.

**2. How much does the T-end rule matter on this corpus?** An independent tangent
(Lepeschkin) T-end is measured on the same median beats and compared with the
study's own annotation. The answer is a small median offset with a
**rate-dependent** component, which therefore does not cancel in the ΔΔ cascade.
"""
from __future__ import annotations

import argparse, sys, warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from dish2ecg.viz import read_wfdb, read_wfdb_annotations   # noqa: E402

#: QT adaptation after a rate step, seconds. Textbook physiology, not measured here.
QT_ADAPTATION_S = (60.0, 120.0)
WINDOW_S = 10.0


def detect_r_peaks(signal, fs):
    """Pan-Tompkins-style R detection. Validated against the table RR, not published."""
    x = np.asarray(signal, float) - np.median(signal)
    ma = lambda v, n: np.convolve(v, np.ones(n) / n, mode="same")
    bp = ma(x, 5) - ma(x, 30)
    integ = ma(np.diff(bp, prepend=bp[0]) ** 2, int(0.08 * fs))
    above = integ > 0.35 * np.percentile(integ, 99)
    starts = np.where(above[1:] & ~above[:-1])[0]
    peaks, refractory = [], int(0.22 * fs)
    for i in starts:
        j1 = min(len(x), i + int(0.15 * fs))
        p = i + int(np.argmax(np.abs(bp[i:j1])))
        if not peaks or p - peaks[-1] > refractory:
            peaks.append(p)
    return np.array(peaks)


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
    return dict(P_onset=on[0], QRS_onset=on[1], J=off[0], T_peak=tpk[0], T_end=off[1])


def main(out: Path, wf: Path, index: Path) -> None:
    warnings.filterwarnings("ignore")
    out.mkdir(parents=True, exist_ok=True)
    sel = pd.read_csv(index)
    fri = lambda qt, rr: qt / (rr / 1000.0) ** (1 / 3)

    rr_rows, dl_rows, failed = [], [], {"raw": 0, "median": 0}
    for _, r in sel.iterrows():
        uid = r["EGREFID"]
        # --- 1. beat-to-beat RR from the raw strip
        if (wf / "raw" / f"{uid}.dat").exists():
            try:
                rec = read_wfdb(wf / "raw" / uid)
                pk = detect_r_peaks(rec.lead("II"), rec.fs)
                rr = np.diff(pk) * 1000.0 / rec.fs
                rr = rr[(rr > 300) & (rr < 2200)]
                if len(rr) >= 3:
                    rr_rows.append(dict(
                        EGREFID=uid, RANDID=r["RANDID"], arm=r["EXTRT"],
                        table_RR=r["RR"], table_QT=r["QT"], n_beats=len(rr) + 1,
                        rr_mean=rr.mean(), rr_sd=rr.std(ddof=1), rr_min=rr.min(),
                        rr_max=rr.max(), rr_range=float(np.ptp(rr)),   # ndarray.ptp() was removed in NumPy 2.0
                        rr_first=rr[0], rr_last=rr[-1]))
            except Exception:
                failed["raw"] += 1
        # --- 2. delineation on the median beat
        if (wf / "medians" / f"{uid}.atr").exists():
            try:
                rec = read_wfdb(wf / "medians" / uid)
                f = fiducials(read_wfdb_annotations(wf / "medians" / f"{uid}.atr"), rec.fs)
                if f is None:
                    continue
                try:
                    y = rec.lead("VCGMAG")
                except KeyError:
                    y = np.abs(rec.lead("II"))
                te = tangent_tend(y, f["T_peak"], rec.fs)
                if np.isfinite(te):
                    dl_rows.append(dict(
                        EGREFID=uid, RANDID=r["RANDID"], arm=r["EXTRT"],
                        table_QT=r["QT"], table_RR=r["RR"], t_peak=f["T_peak"],
                        qt_study=f["T_end"] - f["QRS_onset"],
                        qt_tangent=te - f["QRS_onset"]))
            except Exception:
                failed["median"] += 1

    if not rr_rows:
        raise SystemExit("no raw records parsed - check --wf points at the download")
    print(f"parsed {len(rr_rows)} raw, {len(dl_rows)} median; failures {failed}")
    A = pd.DataFrame(rr_rows)
    A["err_vs_table"] = A.rr_mean - A.table_RR
    A["qtcf_spread"] = fri(A.table_QT, A.rr_min) - fri(A.table_QT, A.rr_max)
    A["drift"] = A.rr_last - A.rr_first
    A.to_csv(out / "waveform_rr_window.csv", index=False)

    B = pd.DataFrame(dl_rows)
    B["dQT"] = B.qt_tangent - B.qt_study
    B["plausible"] = B.qt_tangent.between(250, 650)
    ok = B[B.plausible]
    B.to_csv(out / "waveform_delineation.csv", index=False)

    lo, hi = WINDOW_S / QT_ADAPTATION_S[1], WINDOW_S / QT_ADAPTATION_S[0]
    sl, ic, _, p, _ = stats.linregress(ok.table_RR, ok.dQT)
    print(f"R-detector vs table RR: median {A.err_vs_table.median():+.2f} ms, "
          f"SD {A.err_vs_table.std():.1f}, |err|>20 ms in {100*(A.err_vs_table.abs()>20).mean():.1f}%")
    print(f"window covers {lo:.0%}-{hi:.0%} of one QT adaptation time constant")
    print(f"within-window RR: median SD {A.rr_sd.median():.1f} ms, median range {A.rr_range.median():.1f} ms")
    print(f"implied QTcF spread: median {A.qtcf_spread.median():.1f} ms, "
          f">10 ms in {100*(A.qtcf_spread>10).mean():.1f}% of {len(A)} records")
    print(f"delineation (n={len(ok)}, {100*(~B.plausible).mean():.1f}% of tangent fits rejected): "
          f"median {ok.dQT.median():+.2f} ms, IQR [{ok.dQT.quantile(.25):+.1f}, {ok.dQT.quantile(.75):+.1f}]")
    print(f"  rate dependence: {sl:+.4f} ms per ms RR = {sl*300:+.2f} ms over 300 ms (p={p:.4f})")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=Path("results"))
    ap.add_argument("--wf", type=Path, default=Path("../wfs"))
    ap.add_argument("--index", type=Path, default=Path("../wfs/selected.csv"))
    a = ap.parse_args()
    main(a.out, a.wf, a.index)
