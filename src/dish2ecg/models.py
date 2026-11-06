"""Concentration-response models.

Clinical: the ICH E14 Q&A (R3) linear mixed-effects concentration-QTc model,

    ddQTcF_ij = (intercept + b_i) + slope * C_ij + eps_ij

with a subject random intercept (and optionally a random slope). The primary
readouts are the slope, the predicted effect at a nominated concentration (usually
the geometric mean Cmax), and the two-sided 90% confidence interval of that
prediction - the regulatory decision is whether the upper bound of that interval
exceeds 10 ms.

In vitro: the same model with the random-effect structure swapped to plate / well /
site / platform / cell line, and concentration expressed as multiples of the free
clinical Cmax.
"""
from __future__ import annotations

from dataclasses import dataclass, asdict
from typing import Sequence

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy import stats

from .provenance import Endpoint, ProvenanceError

__all__ = ["ConcentrationResponseFit", "fit_concentration_response", "emax_fit"]

_Z90 = stats.norm.ppf(0.95)  # one-sided 95% = two-sided 90%


def _q90(n_groups: int) -> float:
    """Two-sided 90% quantile.

    A t quantile on ``n_groups - 1`` degrees of freedom, which is the convention in
    regulatory concentration-QTc analysis and matters here: with ~20 subjects the
    normal quantile is anti-conservative. Measured on 200 simulated crossover
    studies, the normal quantile gave 84.5% coverage of the true effect at Cmax
    against a nominal 90%; the t quantile recovers it.
    """
    df = max(int(n_groups) - 1, 1)
    return float(stats.t.ppf(0.95, df))


@dataclass
class ConcentrationResponseFit:
    """Result of a concentration-response fit."""

    endpoint: str
    domain: str
    units: str
    n_obs: int
    n_groups: int
    slope: float
    slope_se: float
    slope_ci90: tuple[float, float]
    intercept: float
    intercept_se: float
    at_concentration: float | None
    effect_at_concentration: float | None
    effect_ci90: tuple[float, float] | None
    random_effects: str
    converged: bool
    boundary_fallback: bool
    correction_cascade: str
    model_formula: str

    def to_row(self) -> dict:
        d = asdict(self)
        for k in ("slope_ci90", "effect_ci90"):
            v = d.pop(k)
            d[k + "_lower"] = None if v is None else float(v[0])
            d[k + "_upper"] = None if v is None else float(v[1])
        return d

    def exceeds_threshold(self, threshold: float = 10.0) -> bool | None:
        """Does the upper 90% bound of the effect at the nominated concentration
        exceed ``threshold``? ``None`` if no concentration was nominated."""
        if self.effect_ci90 is None:
            return None
        return bool(self.effect_ci90[1] > threshold)


def fit_concentration_response(
    df: pd.DataFrame,
    endpoint: Endpoint,
    *,
    concentration_col: str,
    group_col: str,
    random_slope: bool = False,
    at_concentration: float | None = None,
    extra_fixed: Sequence[str] = (),
    require_steps: Sequence[str] = ("rate_correction", "baseline_subtraction",
                                    "control_subtraction"),
    center_at_target: bool = True,
    center_concentration: bool = False,
) -> ConcentrationResponseFit:
    """Linear mixed-effects concentration-response fit.

    Parameters
    ----------
    df
        Records aligned with ``endpoint.values``.
    endpoint
        The double-corrected endpoint (ddQTcF or ddFPDc). Its provenance is checked
        against ``require_steps`` - pass ``require_steps=()`` to fit an endpoint that
        has deliberately had fewer corrections.
    group_col
        The random-effect grouping: subject in the clinic, well or plate in the dish.
    at_concentration
        Concentration at which to predict the effect, in the units of
        ``concentration_col`` (typically the geometric mean Cmax, or 1.0 when the
        column is already in multiples of free Cmax).
    center_at_target
        Fit the model on ``concentration - at_concentration``. The slope is
        unchanged and the predicted effect becomes the intercept, so its standard
        error is read directly instead of through a delta-method combination of
        three covariance entries. This matters: when the random-effects covariance
        is near-singular - common with ~20 subjects - the off-diagonal covariance
        terms are unreliable and the delta-method interval can explode by many
        orders of magnitude while the point estimate stays correct. Leave this on
        unless you need the uncentred parameterisation.
    """
    if require_steps:
        endpoint.require_steps(*require_steps)
    work = df.copy()
    work["_y"] = np.asarray(endpoint.values, dtype=float)
    work["_c"] = pd.to_numeric(work[concentration_col], errors="coerce")
    cols = ["_y", "_c", group_col] + list(extra_fixed)
    work = work.loc[:, cols].dropna()
    if work.empty:
        raise ValueError("no complete records after dropping missing y / concentration")

    if center_at_target and at_concentration is not None:
        cmean = float(at_concentration)
        centering = "at_target"
    elif center_concentration:
        cmean = float(work["_c"].mean())
        centering = "at_mean"
    else:
        cmean = 0.0
        centering = "none"
    work["_c"] = work["_c"] - cmean

    # Internal scaling. Concentrations across a panel of drugs span several orders
    # of magnitude (1.5 ng/mL for dofetilide, 9400 ng/mL for moxifloxacin), which
    # leaves the mixed-model design matrix badly conditioned and makes the
    # optimiser return an unreliable intercept. Fit on a unit-scaled predictor and
    # map the slope back; the prediction and its interval are unaffected.
    cscale = float(work["_c"].std(ddof=0))
    if not np.isfinite(cscale) or cscale <= 0:
        cscale = 1.0
    work["_c"] = work["_c"] / cscale

    rhs = " + ".join(["_c", *extra_fixed]) if extra_fixed else "_c"
    formula = f"_y ~ {rhs}"
    re_formula = "~_c" if random_slope else None
    model = smf.mixedlm(formula, work, groups=work[group_col], re_formula=re_formula)
    try:
        res = model.fit(reml=True, method="lbfgs")
        converged = bool(getattr(res, "converged", True))
    except Exception:
        res = model.fit(reml=True)
        converged = bool(getattr(res, "converged", False))

    # Boundary case. When the between-group variance is estimated at exactly zero -
    # routine for a double-corrected endpoint with an individual rate correction,
    # which removes most of the between-subject variation - statsmodels returns a
    # degenerate intercept with an astronomically large standard error while the
    # slope stays correct. The mixed model has then collapsed to ordinary least
    # squares, so fit that explicitly and report cluster-robust standard errors.
    boundary = False
    try:
        group_var = float(np.max(np.asarray(res.cov_re)))
    except Exception:
        group_var = np.nan
    bad_intercept = (not np.isfinite(res.bse["Intercept"])
                     or float(res.bse["Intercept"]) > 1e4 * max(
                         float(work["_y"].std(ddof=0)), 1e-12))
    if (np.isfinite(group_var) and group_var <= 0) or bad_intercept:
        boundary = True
        ols = smf.ols(formula, work).fit(
            cov_type="cluster", cov_kwds={"groups": work[group_col]})
        res = ols
        converged = True

    q = _q90(int(work[group_col].nunique()))
    slope = float(res.params["_c"]) / cscale
    slope_se = float(res.bse["_c"]) / cscale
    inter = float(res.params["Intercept"])
    inter_se = float(res.bse["Intercept"])

    eff = eff_ci = None
    if at_concentration is not None:
        cc = float(at_concentration) - cmean
        eff = inter + slope * cc
        if abs(cc) < 1e-12:
            # centred at the target: the prediction IS the intercept
            se = inter_se
        else:
            cov = res.cov_params()
            v = (cov.loc["Intercept", "Intercept"]
                 + cc ** 2 * cov.loc["_c", "_c"]
                 + 2 * cc * cov.loc["Intercept", "_c"])
            se = float(np.sqrt(max(v, 0.0)))
        eff_ci = (eff - q * se, eff + q * se)

    return ConcentrationResponseFit(
        endpoint=endpoint.name,
        domain=endpoint.domain,
        units=endpoint.units,
        n_obs=int(len(work)),
        n_groups=int(work[group_col].nunique()),
        slope=slope,
        slope_se=slope_se,
        slope_ci90=(slope - q * slope_se, slope + q * slope_se),
        intercept=inter,
        intercept_se=inter_se,
        at_concentration=None if at_concentration is None else float(at_concentration),
        effect_at_concentration=None if eff is None else float(eff),
        effect_ci90=None if eff_ci is None else (float(eff_ci[0]), float(eff_ci[1])),
        random_effects=(
            f"{group_col}: between-group variance estimated at 0, refitted as OLS "
            "with cluster-robust SE" if boundary else
            f"{group_col}" + (" (intercept+slope)" if random_slope else " (intercept)")),
        converged=converged,
        boundary_fallback=boundary,
        correction_cascade=" | ".join(str(r) for r in endpoint.provenance),
        model_formula=(formula + (f" | random {re_formula}" if re_formula else "")
                       + f" | centering={centering} | predictor scaled by "
                         f"{cscale:.6g} internally"),
    )


def emax_fit(concentration, effect, *, e0: float | None = None) -> dict:
    """Simple Emax fit ``E = E0 + Emax * C / (EC50 + C)`` by least squares.

    Offered for the in-vitro path where saturation within the tested concentration
    range is common. Returns parameter estimates and their standard errors, or
    ``{"converged": False}`` if the optimiser fails.
    """
    from scipy.optimize import curve_fit

    c = np.asarray(concentration, dtype=float)
    y = np.asarray(effect, dtype=float)
    ok = np.isfinite(c) & np.isfinite(y)
    c, y = c[ok], y[ok]
    if c.size < 4:
        return {"converged": False, "reason": "fewer than 4 usable points"}

    if e0 is None:
        def f(x, E0, Emax, EC50):
            return E0 + Emax * x / (EC50 + x)
        p0 = [float(y.min()), float(y.max() - y.min()) or 1.0,
              float(np.median(c[c > 0])) if np.any(c > 0) else 1.0]
        names = ["E0", "Emax", "EC50"]
    else:
        def f(x, Emax, EC50):
            return e0 + Emax * x / (EC50 + x)
        p0 = [float(y.max() - e0) or 1.0,
              float(np.median(c[c > 0])) if np.any(c > 0) else 1.0]
        names = ["Emax", "EC50"]

    try:
        popt, pcov = curve_fit(f, c, y, p0=p0, maxfev=20000)
    except Exception as exc:  # pragma: no cover - optimiser failure path
        return {"converged": False, "reason": str(exc)}
    se = np.sqrt(np.diag(pcov))
    out = {"converged": True, "n": int(c.size)}
    for n, v, s in zip(names, popt, se):
        out[n] = float(v)
        out[n + "_se"] = float(s)
    if e0 is not None:
        out["E0"] = float(e0)
    return out
