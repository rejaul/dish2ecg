"""Heart-rate / beat-rate correction of repolarization intervals.

The same algebra serves both domains:

* clinical  - interval = QT (or J-Tpeak), cycle length = RR
* in vitro  - interval = FPD (or APD90), cycle length = the interbeat interval

Every public function returns an :class:`~dish2ecg.provenance.Endpoint` carrying the
correction record, so the cascade can never be silently lost.

All intervals and cycle lengths are in **milliseconds** on input and output. The
classical formulas are defined with RR in seconds; the conversion is done internally
and is stated in the recorded formula string.
"""
from __future__ import annotations

from typing import Literal

import numpy as np
import pandas as pd

from .provenance import CorrectionRecord, Domain, Endpoint

__all__ = [
    "power_correction",
    "bazett",
    "fridericia",
    "framingham",
    "hodges",
    "jtpeakc_hnatkova2019",
    "jtpeakc_power",
    "fit_power_exponent",
    "fit_individual_exponents",
    "individual_correction",
    "rate_correct",
]

_MS_PER_S = 1000.0


def _as_series(x, name: str) -> pd.Series:
    if isinstance(x, pd.Series):
        return x.astype(float)
    return pd.Series(np.asarray(x, dtype=float), name=name)


def _wrap(values, name, units, domain, record, index) -> Endpoint:
    s = pd.Series(np.asarray(values, dtype=float), index=index, name=name)
    return Endpoint(values=s, name=name, units=units, domain=domain, provenance=(record,))


# ---------------------------------------------------------------------------
# power-law family:  Ic = I / (CL_seconds ** alpha)
# ---------------------------------------------------------------------------
def power_correction(
    interval_ms,
    cycle_length_ms,
    alpha: float,
    *,
    name: str = "Ic",
    domain: Domain = "clinical",
    method: str = "power",
    citation: str | None = None,
    note: str | None = None,
) -> Endpoint:
    """Generic power-law rate correction ``Ic = I / (CL/1000) ** alpha``.

    ``alpha`` = 1/2 is Bazett, ``alpha`` = 1/3 is Fridericia. For in-vitro data the
    same two exponents are conventionally applied to FPD against the interbeat
    interval and are then written FPDcB and FPDcF.
    """
    i = _as_series(interval_ms, "interval")
    cl = _as_series(cycle_length_ms, "cycle_length")
    if len(i) != len(cl):
        raise ValueError(f"length mismatch: interval {len(i)} vs cycle length {len(cl)}")
    with np.errstate(divide="ignore", invalid="ignore"):
        out = i.to_numpy() / np.power(cl.to_numpy() / _MS_PER_S, alpha)
    out = np.where(np.asarray(cl) > 0, out, np.nan)
    rec = CorrectionRecord(
        step="rate_correction",
        method=method,
        formula=f"Ic = I / (CL_ms/1000)**{alpha:g}",
        parameters={"alpha": float(alpha), "interval_units": "ms", "cycle_length_units": "ms"},
        citation=citation,
        note=note,
    )
    idx = i.index if isinstance(interval_ms, pd.Series) else None
    return _wrap(out, name, "ms", domain, rec, idx)


def bazett(interval_ms, cycle_length_ms, *, name="QTcB", domain: Domain = "clinical") -> Endpoint:
    """Bazett correction, ``alpha = 1/2``. Over-corrects at high heart rates."""
    return power_correction(
        interval_ms, cycle_length_ms, 0.5, name=name, domain=domain, method="bazett",
        citation="Bazett HC, Heart 1920;7:353-370",
        note="Over-corrects above ~70 bpm; retained for comparability with older reports.",
    )


def fridericia(interval_ms, cycle_length_ms, *, name="QTcF", domain: Domain = "clinical") -> Endpoint:
    """Fridericia correction, ``alpha = 1/3``. The ICH E14 default."""
    return power_correction(
        interval_ms, cycle_length_ms, 1.0 / 3.0, name=name, domain=domain, method="fridericia",
        citation="Fridericia LS, Acta Med Scand 1920;53:469-486; ICH E14 Q&A",
    )


def framingham(interval_ms, cycle_length_ms, *, name="QTcFram", domain: Domain = "clinical") -> Endpoint:
    """Linear Framingham correction ``QTc = QT + 154 * (1 - RR_s)`` (ms)."""
    i = _as_series(interval_ms, "interval")
    cl = _as_series(cycle_length_ms, "cycle_length")
    out = i.to_numpy() + 154.0 * (1.0 - cl.to_numpy() / _MS_PER_S)
    rec = CorrectionRecord(
        step="rate_correction", method="framingham",
        formula="QTc_ms = QT_ms + 154 * (1 - RR_ms/1000)",
        parameters={"slope_ms_per_s": 154.0},
        citation="Sagie A et al., Am J Cardiol 1992;70:797-801",
    )
    idx = i.index if isinstance(interval_ms, pd.Series) else None
    return _wrap(out, name, "ms", domain, rec, idx)


def hodges(interval_ms, cycle_length_ms, *, name="QTcH", domain: Domain = "clinical") -> Endpoint:
    """Hodges correction ``QTc = QT + 1.75 * (HR - 60)`` (ms)."""
    i = _as_series(interval_ms, "interval")
    cl = _as_series(cycle_length_ms, "cycle_length")
    hr = 60.0 * _MS_PER_S / cl.to_numpy()
    out = i.to_numpy() + 1.75 * (hr - 60.0)
    rec = CorrectionRecord(
        step="rate_correction", method="hodges",
        formula="QTc_ms = QT_ms + 1.75 * (60000/RR_ms - 60)",
        parameters={"slope_ms_per_bpm": 1.75},
        citation="Hodges M et al., J Electrocardiol 1983;16:17-24",
    )
    idx = i.index if isinstance(interval_ms, pd.Series) else None
    return _wrap(out, name, "ms", domain, rec, idx)


def jtpeakc_hnatkova2019(jtpeak_ms, rr_ms, *, name="JTpeakc", domain: Domain = "clinical") -> Endpoint:
    """J-Tpeak rate correction of Vicente et al. 2019.

    ``JTpc = JTp + 0.150 * (1 - RR)`` with every quantity in **seconds**; this
    function takes and returns milliseconds, so the applied constant is 150 ms/s.

    The published derivation uses hysteresis-corrected RR intervals. If the RR you
    pass is the instantaneous RR of the analysed beat, the correction is an
    approximation - see :func:`dish2ecg.cascade.hysteresis_corrected_rr`.
    """
    i = _as_series(jtpeak_ms, "jtpeak")
    cl = _as_series(rr_ms, "rr")
    out = i.to_numpy() + 150.0 * (1.0 - cl.to_numpy() / _MS_PER_S)
    rec = CorrectionRecord(
        step="rate_correction", method="hnatkova2019_jtpeak",
        formula="JTpc_s = JTp_s + 0.150 * (1 - RR_s)  [applied in ms as +150*(1-RR_ms/1000)]",
        parameters={"slope_ms_per_s": 150.0},
        citation="Hnatkova K, Vicente J, Johannesen L, Garnett C, Strauss DG, Stockbridge N, Malik M, Sci Rep 2019;9:15060. doi:10.1038/s41598-019-51491-4",
        note="Derived with hysteresis-corrected RR in 523 healthy subjects.",
    )
    idx = i.index if isinstance(jtpeak_ms, pd.Series) else None
    return _wrap(out, name, "ms", domain, rec, idx)


#: Exponent of the power-law J-Tpeak correction used in the FDA CiPA ECG work.
#: PUBLISHED, not discovered here: Vicente 2019 (Clin Pharmacol Ther 105:943-953,
#: doi:10.1002/cpt.1303) states the prespecified population formula
#: J-Tpeakc = J-Tpeak / (RR/1000 ms)**0.58 and attributes it to Johannesen 2014
#: (Clin Pharmacol Ther 95:501-508, doi:10.1038/clpt.2013.238; not read here).
#: The exponent was read from the PMC XML of Vicente 2019, because the plain-text
#: rendering drops superscripts. Independent check: regressing log(JTP/JTPC) on
#: log(RR) over the 5740 complete records of the study's own derived values in
#: PhysioNet ECGCIPA recovers 0.5799 (pure power law, intercept 0.000, residual
#: SD 0.15 ms).
JTPEAK_POWER_ALPHA_CIPA = 0.58


def jtpeakc_power(
    jtpeak_ms, rr_ms, *, alpha: float = JTPEAK_POWER_ALPHA_CIPA,
    name: str = "JTpeakc", domain: Domain = "clinical",
) -> Endpoint:
    """Power-law J-Tpeak correction ``JTpc = JTp / (RR_s)**alpha``.

    This is the correction used in the FDA's CiPA ECG biomarker work, and it is
    **not** the same as :func:`jtpeakc_hnatkova2019`, which is linear. Over the
    5740 complete records of PhysioNet ECGCIPA the two differ by a median of
    0.44 ms with a standard deviation of 3.56 ms and a maximum of 21.1 ms -
    material against a 10 ms threshold, so the choice has to be deliberate.

    Use this one to reproduce a CiPA-style analysis; use the Hnatkova 2019 linear
    form when following that paper, remembering it expects hysteresis-corrected
    RR.
    """
    return power_correction(
        jtpeak_ms, rr_ms, alpha, name=name, domain=domain, method="jtpeak_power",
        citation=("Vicente J et al., Clin Pharmacol Ther 2019;105:943-953, "
                  "doi:10.1002/cpt.1303: prespecified population formula, exponent "
                  "0.58, attributed there to Johannesen L et al., Clin Pharmacol "
                  "Ther 2014;95:501-508, doi:10.1038/clpt.2013.238; reproduced from "
                  "the study's own derived J-Tpeakc in PhysioNet ECGCIPA v1.0.0"),
        note=("alpha recovered as 0.5799 from 5740 records, intercept 0.000, "
              "residual SD 0.15 ms."),
    )


# ---------------------------------------------------------------------------
# data-driven exponents
# ---------------------------------------------------------------------------
def fit_power_exponent(interval_ms, cycle_length_ms) -> tuple[float, dict]:
    """Fit a pooled study-specific exponent by log-log regression.

    ``log(I) = a + alpha * log(CL_s)``. Fit this on **drug-free** records only
    (baseline and/or placebo), which is how study-specific and population
    corrections are defined.

    Returns ``(alpha, diagnostics)``.
    """
    i = _as_series(interval_ms, "interval")
    cl = _as_series(cycle_length_ms, "cycle_length")
    ok = np.isfinite(i) & np.isfinite(cl) & (i > 0) & (cl > 0)
    if ok.sum() < 3:
        raise ValueError(f"need >=3 usable records to fit an exponent, got {int(ok.sum())}")
    x = np.log(cl.to_numpy()[ok.to_numpy()] / _MS_PER_S)
    y = np.log(i.to_numpy()[ok.to_numpy()])
    alpha, intercept = np.polyfit(x, y, 1)
    resid = y - (intercept + alpha * x)
    ss_tot = float(((y - y.mean()) ** 2).sum())
    diag = {
        "n": int(ok.sum()),
        "intercept_log": float(intercept),
        "r_squared": float(1.0 - (resid ** 2).sum() / ss_tot) if ss_tot > 0 else np.nan,
        "residual_sd_log": float(resid.std(ddof=2)),
    }
    return float(alpha), diag


def fit_individual_exponents(
    df: pd.DataFrame, *, interval_col: str, cycle_length_col: str, subject_col: str,
    min_records: int = 10,
) -> pd.DataFrame:
    """Fit one exponent per subject (or per well) on drug-free records.

    Subjects with fewer than ``min_records`` usable records get ``alpha = NaN`` and
    must fall back to the pooled exponent; the fallback is recorded explicitly by
    :func:`individual_correction`.
    """
    rows = []
    for key, g in df.groupby(subject_col, sort=True):
        try:
            a, d = fit_power_exponent(g[interval_col], g[cycle_length_col])
        except ValueError:
            a, d = np.nan, {"n": len(g), "intercept_log": np.nan,
                            "r_squared": np.nan, "residual_sd_log": np.nan}
        if d["n"] < min_records:
            a = np.nan
        rows.append({subject_col: key, "alpha": a, **d})
    return pd.DataFrame(rows)


def individual_correction(
    df: pd.DataFrame, *, interval_col: str, cycle_length_col: str, subject_col: str,
    exponents: pd.DataFrame, fallback_alpha: float | None = None,
    name: str = "QTcI", domain: Domain = "clinical",
) -> Endpoint:
    """Apply per-subject exponents fitted by :func:`fit_individual_exponents`."""
    amap = exponents.set_index(subject_col)["alpha"].to_dict()
    alphas = df[subject_col].map(amap).astype(float)
    n_fallback = int(alphas.isna().sum())
    if n_fallback and fallback_alpha is None:
        raise ValueError(
            f"{n_fallback} records have no individual exponent and no fallback_alpha "
            "was given."
        )
    if n_fallback:
        alphas = alphas.fillna(float(fallback_alpha))
    cl_s = df[cycle_length_col].to_numpy(dtype=float) / _MS_PER_S
    out = df[interval_col].to_numpy(dtype=float) / np.power(cl_s, alphas.to_numpy())
    rec = CorrectionRecord(
        step="rate_correction", method="individual_power",
        formula="Ic = I / (CL_ms/1000)**alpha_subject",
        parameters={
            "n_subjects_fitted": int(np.isfinite(list(amap.values())).sum()),
            "n_records_on_fallback": n_fallback,
            "fallback_alpha": fallback_alpha,
            "alpha_median": float(np.nanmedian(list(amap.values()))),
        },
        citation="Malik M et al., Drug Saf 2019;42:401-414, doi:10.1007/s40264-018-0736-1 (Implications of Individual QT/RR Profiles, Part 1)",
        note=(f"{n_fallback} records used the pooled fallback exponent."
              if n_fallback else None),
    )
    return _wrap(out, name, "ms", domain, rec, df.index)


_REGISTRY = {
    "bazett": bazett,
    "fridericia": fridericia,
    "framingham": framingham,
    "hodges": hodges,
    "hnatkova2019_jtpeak": jtpeakc_hnatkova2019,
    # deprecated alias: this correction was published by Hnatkova et al., not Vicente
    "vicente2019_jtpeak": jtpeakc_hnatkova2019,
    "jtpeak_power": jtpeakc_power,
}


def rate_correct(
    interval_ms, cycle_length_ms, method: str = "fridericia", *,
    alpha: float | None = None, name: str | None = None, domain: Domain = "clinical",
) -> Endpoint:
    """Dispatch to a named correction. ``method="power"`` requires ``alpha``."""
    if method == "power":
        if alpha is None:
            raise ValueError("method='power' requires alpha")
        return power_correction(interval_ms, cycle_length_ms, alpha,
                                name=name or "Ic", domain=domain)
    try:
        fn = _REGISTRY[method]
    except KeyError:
        raise ValueError(
            f"unknown correction {method!r}; available: "
            f"{sorted(_REGISTRY) + ['power']}"
        ) from None
    kw = {"domain": domain}
    if name:
        kw["name"] = name
    return fn(interval_ms, cycle_length_ms, **kw)


# Deprecated alias. The linear J-Tpeak correction was published by Hnatkova K,
# Vicente J, Johannesen L, Garnett C, Strauss DG, Stockbridge N, Malik M,
# Sci Rep 2019;9:15060 (doi:10.1038/s41598-019-51491-4). Earlier versions of this
# package named it for the second author, which collided with Vicente 2019 CPT.
jtpeakc_vicente2019 = jtpeakc_hnatkova2019
