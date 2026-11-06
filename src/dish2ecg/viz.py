"""Visualising what the cascade measures.

Three kinds of picture, in increasing order of how much they are *data* rather
than *diagram*:

* :func:`plot_ecg_beat` - a **real** median beat from PhysioNet, annotated with
  the study's own fiducial points. Nothing here is drawn by hand.
* :func:`plot_qt_rr` and :func:`plot_concentration_response` - real measured
  values with fitted curves.
* :func:`plot_field_potential` - a **schematic**. No public dataset in this
  project ships hiPSC-CM waveforms; the CiPA myocyte workbook publishes summary
  durations only. The function labels its own output as a schematic and will not
  pretend otherwise.

The WFDB reader here is deliberately minimal (format 16 only) so that the
package gains no dependency for what is a few lines of struct parsing. For
anything beyond reading a median beat, use the ``wfdb`` package.
"""
from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np

__all__ = ["WfdbRecord", "read_wfdb", "read_wfdb_annotations", "ANNOTATION_CODES",
           "confidence_band", "plot_confidence_band", "variance_components",
           "bland_altman", "plot_bland_altman",
           "plot_ecg_beat", "plot_field_potential", "plot_qt_rr",
           "plot_concentration_response", "synthetic_field_potential"]

#: The WFDB annotation codes this module understands.
ANNOTATION_CODES = {1: "N", 24: "p", 27: "t", 39: "(", 40: ")", 63: "aux"}


@dataclass(frozen=True)
class WfdbRecord:
    """A decoded WFDB record. ``signal`` is (samples, leads) in physical units."""

    signal: np.ndarray
    fs: float
    leads: tuple[str, ...]
    units: tuple[str, ...]
    gain: tuple[float, ...]
    baseline: tuple[int, ...]
    adc_resolution: int

    def lead(self, name: str) -> np.ndarray:
        """One lead by name, case-insensitive."""
        up = [l.upper() for l in self.leads]
        if name.upper() not in up:
            raise KeyError(f"no lead {name!r}; have {self.leads}")
        return self.signal[:, up.index(name.upper())]

    @property
    def time_ms(self) -> np.ndarray:
        return np.arange(self.signal.shape[0]) * 1000.0 / self.fs

    def quantisation_uv(self, lead_index: int = 0) -> float:
        """Volts per least-significant bit, in microvolts - the ADC step size."""
        return 1e3 / self.gain[lead_index]


def read_wfdb(record: str | Path) -> WfdbRecord:
    """Read a WFDB format-16 record. ``record`` is the path without extension.

    Physical value = (raw - baseline) / gain, which is the WFDB convention.
    Raises on any format this reader does not implement, rather than guessing.
    """
    record = Path(record)
    head = record.with_suffix(".hea").read_text().strip().splitlines()
    head = [l for l in head if l.strip() and not l.startswith("#")]
    fields = head[0].split()
    nsig, fs = int(fields[1]), float(fields[2])
    leads, units, gains, bases, fmts, res = [], [], [], [], [], []
    for line in head[1:1 + nsig]:
        p = line.split()
        fmts.append(p[1])
        m = re.match(r"([-\d.eE+]+)(?:\(([-\d]+)\))?(?:/(\S+))?", p[2])
        if m is None:
            raise ValueError(f"cannot parse gain field {p[2]!r}")
        unit = m.group(3) or "mV"
        g = float(m.group(1)) or 200.0
        # Headers may quote gain per uV (PhysioNet ECGCIPA) or per mV (ECGRDVQ).
        # Normalise everything to ADC units per mV so `signal` is always in mV.
        if unit.lower() in ("uv", "\u00b5v"):
            g *= 1000.0
            unit = "mV"
        gains.append(g)
        bases.append(int(m.group(2) or 0))
        units.append(unit)
        r_ = int(p[3]) if len(p) > 3 and p[3].lstrip("-").isdigit() else 0
        res.append(r_ or int(p[1]) if p[1] in ("16", "32") else (r_ or 16))
        leads.append(p[-1])
    DTYPE = {"16": "<i2", "32": "<i4"}
    bad = {f for f in fmts if f not in DTYPE}
    if bad:
        raise NotImplementedError(
            f"this minimal reader handles WFDB formats {sorted(DTYPE)} only, found "
            f"{sorted(bad)}; use the 'wfdb' package for other formats")
    if len(set(fmts)) > 1:
        raise NotImplementedError(f"mixed sample formats in one record: {sorted(set(fmts))}")
    raw = np.fromfile(record.with_suffix(".dat"), dtype=DTYPE[fmts[0]])
    raw = raw.reshape(-1, nsig)
    phys = (raw.astype(float) - np.asarray(bases)) / np.asarray(gains)
    return WfdbRecord(signal=phys, fs=fs, leads=tuple(leads), units=tuple(units),
                      gain=tuple(gains), baseline=tuple(bases),
                      adc_resolution=int(res[0]))


def read_wfdb_annotations(path: str | Path) -> list[tuple[int, str]]:
    """Read a WFDB annotation file into ``[(sample, code), ...]``.

    Each 16-bit word carries a 6-bit type code and a 10-bit interval since the
    previous annotation. AUX records (code 63) carry a byte count and a payload,
    which this reader skips rather than decodes.
    """
    words = np.fromfile(Path(path), dtype="<u2")
    out: list[tuple[int, str]] = []
    t, i = 0, 0
    while i < len(words):
        w = int(words[i])
        code, interval = (w >> 10) & 0x3F, w & 0x3FF
        if code == 0 and interval == 0:
            break
        if code == 63:                      # AUX: skip the payload words
            i += 1 + (interval + 1) // 2
            continue
        t += interval
        out.append((t, ANNOTATION_CODES.get(code, str(code))))
        i += 1
    return out


def synthetic_field_potential(duration_ms: float = 352.0, n: int = 4000,
                              span_ms: float = 580.0) -> tuple[np.ndarray, np.ndarray]:
    """A schematic MEA field potential. **Not data** - for explanatory figures.

    Returns ``(time_ms, millivolts)`` with the repolarisation peak at
    ``duration_ms`` after the depolarisation spike.
    """
    t = np.linspace(-60.0, span_ms, n)
    g = lambda c, w, a: a * np.exp(-0.5 * ((t - c) / w) ** 2)
    fp = (g(2.0, 1.6, -1.00) + g(7.5, 4.5, 0.34) + g(95.0, 70.0, -0.045)
          + g(duration_ms, 46.0, 0.26))
    fp[t < -5] = 0.0
    return t, fp


# ---------------------------------------------------------------------------
# plotting
# ---------------------------------------------------------------------------
_FOCAL, _ACC, _WARN, _GREY = "#1b4f9c", "#b3541e", "#9c1b32", "#8a8a8a"


def _axes(ax):
    import matplotlib.pyplot as plt
    if ax is None:
        _, ax = plt.subplots(figsize=(7.2, 3.2))
    return ax


def _span(ax, x0, x1, y, label, color, fontsize=7.5):
    ax.annotate("", xy=(x0, y), xytext=(x1, y),
                arrowprops=dict(arrowstyle="<->", color=color, lw=1.1))
    ax.text((x0 + x1) / 2, y, label, color=color, ha="center", va="bottom",
            fontsize=fontsize, bbox=dict(fc="white", ec="none", pad=0.6))


def plot_ecg_beat(record, *, lead="II", fiducials=None, intervals=None,
                  ax=None, title=None, window_ms=None,
                  show_intervals=("QRS", "JTPEAK", "TPEAKTEND", "QT")):
    """Annotate a median beat with the intervals the cascade consumes.

    ``record``     a :class:`WfdbRecord`, or ``(time_ms, millivolts)``.
    ``fiducials``  ``{"P_onset": ms, "QRS_onset": ms, "J": ms, "T_peak": ms, "T_end": ms}``.
    ``intervals``  measured values in ms, e.g. from the PhysioNet clinical table.
                   Drawn as given; this function never re-measures them.

    Returns the axes. The caller owns the figure.
    """
    ax = _axes(ax)
    if isinstance(record, WfdbRecord):
        t, y = record.time_ms, record.lead(lead)
    else:
        t, y = np.asarray(record[0], float), np.asarray(record[1], float)
    if window_ms is not None:
        m = (t >= window_ms[0]) & (t <= window_ms[1])
        t, y = t[m], y[m]
    ax.plot(t, y, color=_FOCAL, lw=1.25)
    ax.axhline(0.0, color=_GREY, lw=0.6, ls=(0, (4, 3)))
    f = dict(fiducials or {})
    for key, lbl in (("P_onset", "P onset"), ("QRS_onset", "QRS onset"),
                     ("J", "J point"), ("T_peak", "T peak"), ("T_end", "T end")):
        if key in f:
            ax.axvline(f[key], color=_ACC, lw=0.7, ls=(0, (2, 2)), zorder=0)
            ax.text(f[key], float(np.max(y)), lbl, rotation=90, fontsize=6.0,
                    color=_ACC, ha="right", va="top")
    iv = dict(intervals or {})
    lo, hi = float(np.min(y)), float(np.max(y))
    step = 0.09 * (hi - lo)
    base = hi + 0.10 * (hi - lo)
    rows = [("PR", "P_onset", "QRS_onset"), ("QRS", "QRS_onset", "J"),
            ("TPEAKTEND", "T_peak", "T_end"), ("JTPEAK", "J", "T_peak"),
            ("QT", "QRS_onset", "T_end")]
    drawn = [(n, a, b) for n, a, b in rows
             if n in show_intervals and a in f and b in f]
    drawn.sort(key=lambda r: f[r[2]] - f[r[1]])      # shortest span lowest
    step = 0.13 * (hi - lo)
    for k, (name, a, b) in enumerate(drawn):
        lbl = f"{name} {iv[name]:.0f} ms" if name in iv else name
        _span(ax, f[a], f[b], base + k * step, lbl, _ACC)
    k = len(drawn)
    ax.set_xlabel("time (ms)")
    ax.set_ylabel(f"lead {lead} (mV)" if isinstance(record, WfdbRecord) else "mV")
    ax.set_ylim(lo - 0.12 * (hi - lo), base + (k + 0.6) * step)
    if title:
        ax.set_title(title, fontsize=9)
    return ax


def plot_field_potential(fp=None, *, duration_ms=352.0, ax=None, title=None,
                         annotate=True):
    """Annotate a field potential with FPD. Schematic unless ``fp`` is supplied.

    ``fp`` may be ``(time_ms, millivolts)`` from a real recording. When it is
    omitted a schematic is drawn and the axes are labelled as such, because no
    public dataset in this project ships hiPSC-CM waveforms.
    """
    ax = _axes(ax)
    schematic = fp is None
    t, y = synthetic_field_potential(duration_ms) if schematic else (
        np.asarray(fp[0], float), np.asarray(fp[1], float))
    ax.plot(t, y, color=_FOCAL if not schematic else _GREY, lw=1.25)
    ax.axhline(0.0, color=_GREY, lw=0.6, ls=(0, (4, 3)))
    if annotate:
        lo, hi = float(np.min(y)), float(np.max(y))
        for xv in (0.0, duration_ms):
            ax.plot([xv, xv], [lo, hi + 0.22 * (hi - lo)], color=_ACC, lw=0.7,
                    ls=(0, (2, 2)), zorder=0)
        _span(ax, 0.0, duration_ms, hi + 0.24 * (hi - lo),
              f"FPD {duration_ms:.0f} ms", _ACC)
        ax.annotate("depolarisation\nspike (fast Na$^+$)", xy=(3, lo * 0.92),
                    xytext=(0.18, 0.12), textcoords="axes fraction", fontsize=6.5,
                    color=_GREY, arrowprops=dict(arrowstyle="-", color=_GREY, lw=0.6))
        ax.annotate("repolarisation-wave\npeak (IKr)", xy=(duration_ms, hi * 0.96),
                    xytext=(0.52, 0.30), textcoords="axes fraction", fontsize=6.5,
                    color=_GREY, ha="center",
                    arrowprops=dict(arrowstyle="-", color=_GREY, lw=0.6))
        ax.set_ylim(lo - 0.12 * (hi - lo), hi + 0.42 * (hi - lo))
    ax.set_xlabel("time from depolarisation (ms)")
    ax.set_ylabel("extracellular potential (mV)")
    ax.set_title(title or ("Field potential (SCHEMATIC - not a recording)"
                           if schematic else "Field potential"), fontsize=9)
    return ax


def plot_qt_rr(rr_ms, qt_ms, *, alphas=(("Bazett", 0.5), ("Fridericia", 1 / 3)),
               fitted_alpha=None, ax=None, title=None):
    """Measured QT against RR, with correction curves through the data centroid.

    Each curve is ``QT = a * RR_s**alpha`` with ``a`` set so the curve passes
    through the geometric centroid of the plotted points, which is what makes
    the exponents visually comparable.
    """
    ax = _axes(ax)
    rr = np.asarray(rr_ms, float); qt = np.asarray(qt_ms, float)
    ok = np.isfinite(rr) & np.isfinite(qt) & (rr > 0) & (qt > 0)
    rr, qt = rr[ok], qt[ok]
    ax.plot(rr, qt, "o", ms=2.6, mfc=_FOCAL, mec="none", alpha=0.38,
            label=f"n = {rr.size}")
    gx, gy = np.exp(np.mean(np.log(rr))), np.exp(np.mean(np.log(qt)))
    grid = np.linspace(rr.min(), rr.max(), 200)
    curves = list(alphas) + ([("fitted", float(fitted_alpha))] if fitted_alpha else [])
    for (name, a), col in zip(curves, (_GREY, _ACC, _WARN)):
        ax.plot(grid, gy * (grid / gx) ** a, color=col, lw=1.3,
                label=f"{name} $\\alpha$={a:.3f}")
    ax.set_xlabel("RR (ms)"); ax.set_ylabel("QT (ms)")
    ax.legend(fontsize=6.5, frameon=False)
    if title:
        ax.set_title(title, fontsize=9)
    return ax


def plot_concentration_response(conc, effect, *, slope=None, intercept=None,
                                at_concentration=None, effect_at=None, ci=None,
                                log_x=False, ax=None, title=None,
                                xlabel=None, ylabel="$\\Delta\\Delta$QTcF (ms)",
                                threshold=10.0):
    """Concentration-response points with the fitted line and the decision point.

    ``threshold`` draws the regulatory line (10 ms for clinical QTc). Pass
    ``None`` to omit it - there is no agreed in-vitro equivalent.
    """
    ax = _axes(ax)
    c = np.asarray(conc, float); e = np.asarray(effect, float)
    ok = np.isfinite(c) & np.isfinite(e)
    ax.plot(c[ok], e[ok], "o", ms=3.0, mfc=_FOCAL, mec="none", alpha=0.40)
    if slope is not None and intercept is not None and at_concentration is not None:
        grid = np.linspace(np.nanmin(c[ok]), np.nanmax(c[ok]), 200)
        ax.plot(grid, intercept + slope * (grid - at_concentration), color=_ACC, lw=1.4)
    if threshold is not None:
        ax.axhline(threshold, color=_WARN, lw=0.9, ls=(0, (3, 2)))
        ax.text(ax.get_xlim()[1], threshold, f" {threshold:.0f} ms", color=_WARN,
                fontsize=6.5, va="bottom", ha="right")
    ax.axhline(0.0, color=_GREY, lw=0.6)
    if at_concentration is not None and effect_at is not None:
        ax.plot([at_concentration], [effect_at], "o", ms=6, mfc="white",
                mec=_ACC, mew=1.6, zorder=5)
        if ci is not None:
            ax.plot([at_concentration, at_concentration], list(ci), color=_ACC, lw=1.6)
            ax.annotate(f"{effect_at:.1f} [{ci[0]:.1f}, {ci[1]:.1f}]",
                        xy=(at_concentration, ci[1]), xytext=(0.04, 0.93),
                        textcoords="axes fraction", color=_ACC, fontsize=6.8,
                        arrowprops=dict(arrowstyle="-", color=_ACC, lw=0.6))
    if log_x:
        ax.set_xscale("log")
    ax.set_xlabel(xlabel or "concentration")
    ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title, fontsize=9)
    return ax


# ---------------------------------------------------------------------------
# Uncertainty: confidence bands and method comparison
# ---------------------------------------------------------------------------

def confidence_band(x, *, slope, intercept, slope_se, intercept_se, centre,
                    df, level=0.90):
    """Pointwise confidence band for a linear fit centred at ``centre``.

    The cascade's concentration-response models are fitted with the predictor
    centred at the prediction point, so the intercept *is* the effect there and
    the slope/intercept covariance is approximately zero at that point. The band
    is therefore

        se(x) = sqrt(intercept_se**2 + ((x - centre) * slope_se)**2)

    with a Student-t quantile on ``df`` degrees of freedom. Pass
    ``df = n_groups - 1`` to reproduce the interval the model layer reports, so
    that a plotted band and a quoted effect can never disagree.

    Returns ``(fit, lower, upper)``.
    """
    import numpy as _np
    from scipy import stats as _stats

    x = _np.asarray(x, dtype=float)
    fit = intercept + slope * (x - centre)
    se = _np.sqrt(intercept_se ** 2 + ((x - centre) * slope_se) ** 2)
    t = _stats.t.ppf(0.5 + level / 2.0, df)
    return fit, fit - t * se, fit + t * se


def plot_confidence_band(x, *, slope, intercept, slope_se, intercept_se, centre,
                         df, level=0.90, ax=None, color="#2B5D9B", label=None,
                         mark_centre=True, centre_label=None):
    """Draw a fitted line with its confidence band, marking the prediction point.

    ``centre`` is normally the geometric-mean observed Cmax; it is marked because
    every effect this project quotes is read off at that point and nowhere else.
    """
    import matplotlib.pyplot as _plt
    import numpy as _np

    ax = ax or _plt.gca()
    x = _np.asarray(x, dtype=float)
    fit, lo, hi = confidence_band(x, slope=slope, intercept=intercept,
                                  slope_se=slope_se, intercept_se=intercept_se,
                                  centre=centre, df=df, level=level)
    ax.fill_between(x, lo, hi, color=color, alpha=0.17, lw=0,
                    label=None if label is None else f"{int(level * 100)} % band")
    ax.plot(x, fit, color=color, lw=1.4, label=label)
    if mark_centre:
        _, c_lo, c_hi = confidence_band([centre], slope=slope, intercept=intercept,
                                        slope_se=slope_se, intercept_se=intercept_se,
                                        centre=centre, df=df, level=level)
        ax.axvline(centre, color="#8A8A8A", lw=0.7, ls=(0, (2, 2)), zorder=0)
        ax.errorbar([centre], [intercept], yerr=[[intercept - c_lo[0]], [c_hi[0] - intercept]],
                    fmt="o", ms=4.5, color=color, capsize=2.4, lw=1.1, zorder=5,
                    mfc="white", mew=1.1)
        if centre_label:
            ax.annotate(centre_label, xy=(centre, intercept), xytext=(6, -2),
                        textcoords="offset points", fontsize=7, color=color, va="top")
    return ax


def variance_components(diff, group):
    """Unbalanced one-way random-effects decomposition of repeated differences.

    Returns ``(var_between, var_within)``. This is the decomposition Bland and
    Altman prescribe when each subject (or drug) contributes several paired
    measurements: using the raw SD of differences instead treats clustered
    observations as independent and understates the limits of agreement.
    """
    import numpy as _np
    import pandas as _pd

    df = _pd.DataFrame({"d": _np.asarray(diff, dtype=float), "g": _np.asarray(group)}).dropna()
    k, n = df.g.nunique(), len(df)
    if k < 2 or n <= k:
        return 0.0, float(df.d.var(ddof=1))
    means = df.groupby("g").d.mean()
    sizes = df.groupby("g").size()
    ms_b = float((sizes * (means - df.d.mean()) ** 2).sum()) / (k - 1)
    ms_w = float(((df.d - df.g.map(means)) ** 2).sum()) / (n - k)
    n0 = (n - (sizes ** 2).sum() / n) / (k - 1)
    return max((ms_b - ms_w) / n0, 0.0), ms_w


def bland_altman(a, b, *, group=None):
    """Bias and limits of agreement for two methods measuring the SAME quantity.

    ``a`` and ``b`` must be the same observable in the same units; the function
    does not and cannot check that. Comparing endpoints that merely resemble one
    another (a field potential duration against a QT interval, say) produces a
    mean axis that estimates nothing - see ``docs/PROJECT_REPORT.md`` section 15.

    With ``group``, the SD of differences is rebuilt from the between- and
    within-cluster variance components.

    Returns a dict with ``mean``, ``diff``, ``bias``, ``sd``, ``loa_lower``,
    ``loa_upper``, ``prop_bias_slope`` and ``n``.
    """
    import numpy as _np

    a = _np.asarray(a, dtype=float)
    b = _np.asarray(b, dtype=float)
    ok = _np.isfinite(a) & _np.isfinite(b)
    a, b = a[ok], b[ok]
    diff, mean = a - b, (a + b) / 2.0
    if group is None:
        sd = float(_np.std(diff, ddof=1))
        var_b, var_w = _np.nan, sd ** 2
    else:
        g = _np.asarray(group)[ok]
        var_b, var_w = variance_components(diff, g)
        sd = float(_np.sqrt(var_b + var_w))
    bias = float(_np.mean(diff))
    slope = float(_np.polyfit(mean, diff, 1)[0]) if len(diff) > 2 else _np.nan
    return dict(mean=mean, diff=diff, bias=bias, sd=sd,
                var_between=var_b, var_within=var_w,
                loa_lower=bias - 1.96 * sd, loa_upper=bias + 1.96 * sd,
                prop_bias_slope=slope, n=int(len(diff)))


def plot_bland_altman(a, b, *, group=None, ax=None, color="#2B5D9B", point_size=7):
    """Bland-Altman plot: difference against mean, with bias and limits drawn."""
    import matplotlib.pyplot as _plt

    ax = ax or _plt.gca()
    s = bland_altman(a, b, group=group)
    ax.scatter(s["mean"], s["diff"], s=point_size, color=color, alpha=0.45, lw=0)
    ax.axhline(s["bias"], color="#1a1a1a", lw=1.0)
    ax.axhline(s["loa_upper"], color="#B3242B", lw=0.9, ls=(0, (4, 2)))
    ax.axhline(s["loa_lower"], color="#B3242B", lw=0.9, ls=(0, (4, 2)))
    ax.axhline(0, color="#8A8A8A", lw=0.6, ls=":")
    ax.set_xlabel("mean of the two methods")
    ax.set_ylabel("difference")
    return ax, s
