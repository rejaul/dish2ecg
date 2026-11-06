"""Putting the dish and the clinic on one scale.

A field potential duration is not a QT interval, and this module never converts one
into the other. What it does is put both on scales that *are* comparable:

* concentration on a shared scale - multiples of the free (unbound) clinical Cmax;
* effect on a shared scale - percent change from the corrected baseline, or a
  standardized effect (change divided by the control-arm standard deviation);
* decision on a shared scale - does the upper 90% bound at 1x free Cmax exceed the
  threshold of regulatory concern?
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from .models import ConcentrationResponseFit
from .provenance import CorrectionRecord, Endpoint

__all__ = ["to_free_cmax_multiples", "to_percent_change", "to_standardized_effect",
           "decision_agreement", "BridgeComparison", "compare_slopes"]


def to_free_cmax_multiples(
    concentration_molar, free_cmax_molar: float, *, name: str = "folds_free_cmax"
) -> pd.Series:
    """Express a concentration as multiples of the free clinical Cmax."""
    if not np.isfinite(free_cmax_molar) or free_cmax_molar <= 0:
        raise ValueError(f"free_cmax_molar must be positive and finite, got {free_cmax_molar!r}")
    c = pd.Series(np.asarray(concentration_molar, dtype=float))
    return (c / float(free_cmax_molar)).rename(name)


def to_percent_change(endpoint: Endpoint, reference_ms, *, name: str | None = None) -> Endpoint:
    """Convert an absolute change (ms) into percent of the corrected reference value.

    ``reference_ms`` is the corrected baseline level of the same records - QTcF at
    baseline in the clinic, FPDc at baseline in the dish.
    """
    if endpoint.units != "ms":
        raise ValueError(f"expected ms, got {endpoint.units!r}")
    ref = pd.Series(np.asarray(reference_ms, dtype=float), index=endpoint.values.index)
    with np.errstate(divide="ignore", invalid="ignore"):
        pct = 100.0 * endpoint.values.to_numpy() / ref.to_numpy()
    pct = pd.Series(np.where(ref.to_numpy() > 0, pct, np.nan), index=endpoint.values.index)
    rec = CorrectionRecord(
        step="scale_transform", method="percent_of_corrected_baseline",
        formula="pct = 100 * delta_ms / corrected_baseline_ms",
        parameters={"reference_median_ms": float(np.nanmedian(ref))},
        note="Percent change is the shared scale for dish-to-clinic comparison; "
             "raw milliseconds are not comparable between a monolayer and a heart.",
    )
    nm = name or (endpoint.name + "_pct")
    return endpoint.derive(pct.rename(nm), rec, name=nm, units="percent")


def to_standardized_effect(endpoint: Endpoint, control_sd: float, *, name=None) -> Endpoint:
    """Divide an effect by the control-arm standard deviation of the same endpoint."""
    if not np.isfinite(control_sd) or control_sd <= 0:
        raise ValueError(f"control_sd must be positive and finite, got {control_sd!r}")
    vals = endpoint.values / float(control_sd)
    rec = CorrectionRecord(
        step="scale_transform", method="standardized_effect",
        formula="z = delta / sd(control arm)",
        parameters={"control_sd": float(control_sd)},
    )
    nm = name or (endpoint.name + "_z")
    return endpoint.derive(vals.rename(nm), rec, name=nm, units="sd_units")


@dataclass
class BridgeComparison:
    drug: str
    invitro_slope: float
    invitro_slope_se: float
    clinical_slope: float
    clinical_slope_se: float
    invitro_effect_at_1x: float | None
    clinical_effect_at_cmax: float | None
    invitro_units: str
    clinical_units: str
    slope_ratio: float | None
    slope_ratio_ci90: tuple[float, float] | None


def compare_slopes(
    drug: str, invitro: ConcentrationResponseFit, clinical: ConcentrationResponseFit
) -> BridgeComparison:
    """Compare an in-vitro and a clinical concentration-response slope.

    The ratio is reported with a 90% interval from the delta method on the log ratio.
    It is a *relative sensitivity*, not a conversion factor, and is only meaningful
    when both slopes were fitted on the same concentration scale (multiples of free
    Cmax) and the same effect scale.
    """
    ratio = ci = None
    if clinical.slope != 0 and np.isfinite(invitro.slope) and np.isfinite(clinical.slope):
        if invitro.slope > 0 and clinical.slope > 0:
            ratio = invitro.slope / clinical.slope
            rel = np.sqrt((invitro.slope_se / invitro.slope) ** 2
                          + (clinical.slope_se / clinical.slope) ** 2)
            from scipy import stats as _st
            z = _st.norm.ppf(0.95)
            ci = (float(ratio * np.exp(-z * rel)), float(ratio * np.exp(z * rel)))
        else:
            ratio = invitro.slope / clinical.slope
    return BridgeComparison(
        drug=drug,
        invitro_slope=invitro.slope, invitro_slope_se=invitro.slope_se,
        clinical_slope=clinical.slope, clinical_slope_se=clinical.slope_se,
        invitro_effect_at_1x=invitro.effect_at_concentration,
        clinical_effect_at_cmax=clinical.effect_at_concentration,
        invitro_units=invitro.units, clinical_units=clinical.units,
        slope_ratio=None if ratio is None else float(ratio),
        slope_ratio_ci90=ci,
    )


def decision_agreement(
    invitro: ConcentrationResponseFit,
    clinical: ConcentrationResponseFit,
    *,
    invitro_threshold: float,
    clinical_threshold: float = 10.0,
) -> dict:
    """Agreement of the two decisions at their own thresholds.

    The clinical threshold is the 10 ms of ICH E14. There is **no** agreed in-vitro
    equivalent, so ``invitro_threshold`` must be supplied explicitly and its
    justification recorded by the caller.
    """
    iv = invitro.exceeds_threshold(invitro_threshold)
    cl = clinical.exceeds_threshold(clinical_threshold)
    return {
        "invitro_positive": iv,
        "clinical_positive": cl,
        "agree": None if (iv is None or cl is None) else bool(iv == cl),
        "invitro_threshold": float(invitro_threshold),
        "invitro_threshold_units": invitro.units,
        "clinical_threshold": float(clinical_threshold),
        "clinical_threshold_units": clinical.units,
        "note": "The in-vitro threshold is a user-supplied convention, not a "
                "regulatory criterion.",
    }
