"""Adapter for microelectrode-array and voltage-dye exports.

Instrument vendors each ship their own wide CSV/XLSX layout. Rather than guess, this
module takes an explicit column mapping and normalises it to the long schema the
package works in. ``MEA_COLUMN_ALIASES`` holds the aliases seen in the wild as a
starting point; anything not matched must be mapped by hand, and the reader says
which columns it could not resolve.

Canonical long schema
---------------------
plate, well, site, platform, cell_line, drug, concentration (mol/L),
replicate, beat_rate_hz or interbeat_interval_ms, fpd_ms and/or apd90_ms,
arrhythmia_event (0/1), event_class, is_baseline, is_vehicle
"""
from __future__ import annotations

from typing import Mapping, Sequence

import numpy as np
import pandas as pd

__all__ = ["CANONICAL_COLUMNS", "MEA_COLUMN_ALIASES", "read_long_mea",
           "interbeat_interval_from_rate"]

CANONICAL_COLUMNS: tuple[str, ...] = (
    "plate", "well", "site", "platform", "cell_line", "drug", "concentration",
    "replicate", "interbeat_interval_ms", "fpd_ms", "apd90_ms",
    "arrhythmia_event", "event_class", "is_baseline", "is_vehicle",
)

#: Lower-cased aliases -> canonical name. Extend freely; nothing is guessed beyond
#: an exact (case- and separator-insensitive) match.
MEA_COLUMN_ALIASES: dict[str, str] = {
    "plateid": "plate", "plate_id": "plate", "plate barcode": "plate",
    "wellid": "well", "well_id": "well", "well name": "well", "wellname": "well",
    "siteid": "site", "site_id": "site", "lab": "site", "laboratory": "site",
    "instrument": "platform", "device": "platform", "system": "platform",
    "cellline": "cell_line", "cell line": "cell_line", "cell_type": "cell_line",
    "celltype": "cell_line", "lot": "cell_line",
    "compound": "drug", "drugname": "drug", "drug_name": "drug", "treatment": "drug",
    "conc": "concentration", "conc_m": "concentration",
    "concentration_m": "concentration", "dose": "concentration",
    "rep": "replicate", "replicate_id": "replicate",
    "ibi": "interbeat_interval_ms", "ibi_ms": "interbeat_interval_ms",
    "interbeat_interval": "interbeat_interval_ms",
    "beat period (ms)": "interbeat_interval_ms",
    "fpd": "fpd_ms", "fpd (ms)": "fpd_ms", "field_potential_duration": "fpd_ms",
    "fieldpotentialduration": "fpd_ms",
    "apd90": "apd90_ms", "apd90 (ms)": "apd90_ms", "apd_90": "apd90_ms",
    "ead": "arrhythmia_event", "arrhythmia": "arrhythmia_event",
    "cellrhythmia": "arrhythmia_event",
    "ead_type": "event_class", "type_of_eads": "event_class",
    "arrhythmia_type": "event_class",
}


def _norm(name: str) -> str:
    return str(name).strip().lower().replace("-", " ").replace("_", " ")


def interbeat_interval_from_rate(beat_rate, *, unit: str = "hz") -> pd.Series:
    """Convert a beat rate to an interbeat interval in milliseconds."""
    r = pd.to_numeric(pd.Series(beat_rate), errors="coerce").astype(float)
    if unit.lower() in ("hz", "beats/s", "bps"):
        per_s = r
    elif unit.lower() in ("bpm", "beats/min", "min-1"):
        per_s = r / 60.0
    else:
        raise ValueError(f"unknown beat-rate unit {unit!r}; use 'hz' or 'bpm'")
    with np.errstate(divide="ignore", invalid="ignore"):
        out = 1000.0 / per_s
    return pd.Series(np.where(per_s > 0, out, np.nan), index=r.index,
                     name="interbeat_interval_ms")


def read_long_mea(
    source, *, mapping: Mapping[str, str] | None = None,
    beat_rate_col: str | None = None, beat_rate_unit: str = "hz",
    sheet_name=0, strict: bool = True,
) -> pd.DataFrame:
    """Read an instrument export into the canonical long schema.

    Parameters
    ----------
    source
        Path to a CSV or Excel file, or a DataFrame.
    mapping
        Explicit ``{source_column: canonical_column}`` overrides, applied after the
        alias table.
    strict
        If True, raise when a required column (``drug``, ``concentration``, and at
        least one of ``fpd_ms``/``apd90_ms``) cannot be resolved.

    Returns
    -------
    DataFrame
        Canonical columns first, unresolved source columns kept with their original
        names. ``df.attrs['unmapped_columns']`` lists what was not recognised.
    """
    if isinstance(source, pd.DataFrame):
        raw = source.copy()
    else:
        s = str(source)
        raw = (pd.read_excel(source, sheet_name=sheet_name)
               if s.lower().endswith((".xlsx", ".xls")) else pd.read_csv(source))

    ren: dict[str, str] = {}
    for c in raw.columns:
        n = _norm(c)
        if n in MEA_COLUMN_ALIASES:
            ren[c] = MEA_COLUMN_ALIASES[n]
        elif n.replace(" ", "_") in CANONICAL_COLUMNS:
            ren[c] = n.replace(" ", "_")
    if mapping:
        ren.update(dict(mapping))
    df = raw.rename(columns=ren)

    if beat_rate_col is not None and "interbeat_interval_ms" not in df.columns:
        src = ren.get(beat_rate_col, beat_rate_col)
        if src not in df.columns:
            raise KeyError(f"beat_rate_col {beat_rate_col!r} not in the table")
        df["interbeat_interval_ms"] = interbeat_interval_from_rate(
            df[src], unit=beat_rate_unit).to_numpy()

    present = [c for c in CANONICAL_COLUMNS if c in df.columns]
    unmapped = [c for c in df.columns if c not in CANONICAL_COLUMNS]

    missing_required = [c for c in ("drug", "concentration") if c not in df.columns]
    if "fpd_ms" not in df.columns and "apd90_ms" not in df.columns:
        missing_required.append("fpd_ms or apd90_ms")
    if missing_required and strict:
        raise KeyError(
            f"could not resolve required column(s) {missing_required}. "
            f"Source columns: {list(raw.columns)}. Pass mapping={{...}} to map them "
            "explicitly; nothing is inferred beyond the alias table."
        )

    out = df.loc[:, present + unmapped]
    out.attrs["unmapped_columns"] = unmapped
    out.attrs["applied_mapping"] = ren
    out.attrs["missing_required"] = missing_required
    return out
