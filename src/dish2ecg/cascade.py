"""The baseline / vehicle / placebo subtraction cascade.

Clinical path   : QTc -> dQTcF (minus time-matched baseline) -> ddQTcF (minus placebo)
In-vitro path   : FPDc -> dFPDc (minus pre-dose baseline)    -> ddFPDc (minus vehicle)

Both run through the same three functions. Only the column names and the identity of
the control arm differ.
"""
from __future__ import annotations

from typing import Hashable, Sequence

import numpy as np
import pandas as pd

from .provenance import CorrectionRecord, Domain, Endpoint, ProvenanceError

__all__ = [
    "hysteresis_corrected_rr",
    "subtract_baseline",
    "subtract_control",
    "double_delta",
]


def hysteresis_corrected_rr(
    rr_history_ms: Sequence[float], tau_s: float = 120.0, dt_s: float | None = None
) -> float:
    """Exponentially weighted RR over a preceding history window.

    QT follows a change in heart rate with a lag of roughly two minutes, so the RR
    that governs the measured QT is not the instantaneous RR of that beat. Supply the
    RR intervals preceding the analysed beat, oldest first.

    Notes
    -----
    This requires beat-to-beat RR history. The PhysioNet ECGRDVQ / ECGDMMLD clinical
    tables publish one summary RR per extracted ECG replicate and therefore do **not**
    support hysteresis correction; the waveform files would have to be reprocessed.
    The package records that fact rather than silently applying an uncorrected RR as
    if it were hysteresis-corrected.
    """
    rr = np.asarray(rr_history_ms, dtype=float)
    rr = rr[np.isfinite(rr)]
    if rr.size == 0:
        return float("nan")
    if dt_s is None:
        dt_s = float(np.nanmean(rr)) / 1000.0
    age_s = np.arange(rr.size - 1, -1, -1, dtype=float) * dt_s
    w = np.exp(-age_s / float(tau_s))
    return float(np.sum(w * rr) / np.sum(w))


def _require(df: pd.DataFrame, cols: Sequence[str]) -> None:
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise KeyError(f"missing column(s) {missing}; have {list(df.columns)}")


def subtract_baseline(
    df: pd.DataFrame,
    endpoint: Endpoint,
    *,
    subject_col: str,
    arm_col: str,
    baseline_mask: pd.Series,
    time_col: str | None = None,
    match_on_time: bool = False,
    name: str | None = None,
) -> Endpoint:
    """Subtract each subject's (or well's) own baseline value within each arm.

    Parameters
    ----------
    baseline_mask
        Boolean series, True for the records that constitute the baseline. For a
        crossover ECG study this is the pre-dose day or pre-dose timepoints of the
        same period; for a plate this is the pre-compound recording of the same well.
    match_on_time
        If True the baseline is matched timepoint by timepoint (requires
        ``time_col``), giving a time-matched baseline. If False the subject-by-arm
        mean of the baseline records is used.

    Returns
    -------
    Endpoint
        ``delta`` values with a ``baseline_subtraction`` record appended.
    """
    _require(df, [subject_col, arm_col] + ([time_col] if time_col else []))
    if len(df) != len(endpoint.values):
        raise ProvenanceError("df and endpoint have different lengths")
    if endpoint.has_step("baseline_subtraction"):
        raise ProvenanceError(f"{endpoint.name!r} already had a baseline subtracted")

    work = df.copy()
    work["_v"] = np.asarray(endpoint.values, dtype=float)
    mask = np.asarray(baseline_mask, dtype=bool)
    base = work.loc[mask]
    if base.empty:
        raise ValueError("baseline_mask selected no records")

    if match_on_time:
        if time_col is None:
            raise ValueError("match_on_time=True requires time_col")
        keys = [subject_col, arm_col, time_col]
        method = "time_matched_baseline"
        formula = "delta = value - baseline(subject, arm, timepoint)"
    else:
        keys = [subject_col, arm_col]
        method = "mean_baseline_within_arm"
        formula = "delta = value - mean(baseline(subject, arm))"

    bmap = base.groupby(keys, sort=False)["_v"].mean()
    idx = pd.MultiIndex.from_frame(work[keys]) if len(keys) > 1 else pd.Index(work[keys[0]])
    bvals = pd.Series(idx.map(bmap), index=work.index, dtype=float)
    n_missing = int(bvals.isna().sum() - np.asarray(endpoint.values.isna()).sum())
    delta = work["_v"] - bvals

    rec = CorrectionRecord(
        step="baseline_subtraction",
        method=method,
        formula=formula,
        parameters={"keys": keys, "n_baseline_records": int(mask.sum()),
                    "n_records_without_baseline": max(n_missing, 0)},
        citation="ICH E14 Q&A (R3), section on baseline definition",
        note=(f"{n_missing} records had no matching baseline and are NaN."
              if n_missing > 0 else None),
    )
    new_name = name or ("d" + endpoint.name)
    return endpoint.derive(delta.rename(new_name), rec, name=new_name)


def subtract_control(
    df: pd.DataFrame,
    endpoint: Endpoint,
    *,
    arm_col: str,
    control_arm: Hashable,
    time_col: str | None = None,
    subject_col: str | None = None,
    paired: bool = False,
    name: str | None = None,
) -> Endpoint:
    """Subtract the control arm: placebo in the clinic, vehicle in the dish.

    Parameters
    ----------
    paired
        If True the control is subtracted within subject (valid for a crossover
        design where the same subject received placebo). If False the control-arm
        mean at the matching timepoint is subtracted.
    """
    _require(df, [arm_col] + ([time_col] if time_col else []) +
                 ([subject_col] if subject_col else []))
    if not endpoint.has_step("baseline_subtraction"):
        raise ProvenanceError(
            f"{endpoint.name!r} has no baseline subtraction; the double-delta endpoint "
            "requires baseline subtraction first. Applied: " + str(endpoint.steps)
        )
    if endpoint.has_step("control_subtraction"):
        raise ProvenanceError(f"{endpoint.name!r} already had a control arm subtracted")

    work = df.copy()
    work["_v"] = np.asarray(endpoint.values, dtype=float)
    ctrl = work.loc[work[arm_col] == control_arm]
    if ctrl.empty:
        raise ValueError(
            f"control arm {control_arm!r} not found in {arm_col!r}; "
            f"present: {sorted(work[arm_col].dropna().unique())}"
        )

    keys: list[str] = []
    if paired:
        if subject_col is None:
            raise ValueError("paired=True requires subject_col")
        keys.append(subject_col)
    if time_col is not None:
        keys.append(time_col)
    if not keys:
        cmap = pd.Series({(): ctrl["_v"].mean()})
        cvals = pd.Series(ctrl["_v"].mean(), index=work.index, dtype=float)
        method = "overall_control_mean"
    else:
        cmap = ctrl.groupby(keys, sort=False)["_v"].mean()
        idx = pd.MultiIndex.from_frame(work[keys]) if len(keys) > 1 else pd.Index(work[keys[0]])
        cvals = pd.Series(idx.map(cmap), index=work.index, dtype=float)
        method = ("paired_time_matched_control" if paired and time_col
                  else "paired_control" if paired else "time_matched_control_mean")

    dd = work["_v"] - cvals
    rec = CorrectionRecord(
        step="control_subtraction",
        method=method,
        formula="dd = delta - delta(control arm" + (", " + ", ".join(keys) if keys else "") + ")",
        parameters={"control_arm": str(control_arm), "keys": keys,
                    "n_control_records": int(len(ctrl))},
        citation="ICH E14 Q&A (R3), placebo-corrected change from baseline",
        note=("Vehicle control in the in-vitro domain."
              if endpoint.domain == "invitro" else None),
    )
    new_name = name or ("d" + endpoint.name)
    out = endpoint.derive(dd.rename(new_name), rec, name=new_name)
    # drop the control arm's own (identically zero on the paired path) records
    return out


def double_delta(
    df: pd.DataFrame,
    endpoint: Endpoint,
    *,
    subject_col: str,
    arm_col: str,
    baseline_mask: pd.Series,
    control_arm: Hashable,
    time_col: str | None = None,
    match_on_time: bool = True,
    paired: bool = True,
    name: str | None = None,
) -> Endpoint:
    """Run baseline then control subtraction in one call.

    Produces ddQTcF in the clinical domain and ddFPDc in the in-vitro domain.
    """
    d = subtract_baseline(
        df, endpoint, subject_col=subject_col, arm_col=arm_col,
        baseline_mask=baseline_mask, time_col=time_col, match_on_time=match_on_time,
    )
    dd = subtract_control(
        df, d, arm_col=arm_col, control_arm=control_arm, time_col=time_col,
        subject_col=subject_col, paired=paired,
        name=name or ("dd" + endpoint.name),
    )
    return dd
