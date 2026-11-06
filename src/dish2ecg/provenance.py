"""Provenance records and the guard rails that travel with every derived endpoint.

The central rule of this package: a derived value is meaningless without the list of
corrections that produced it, and an in-vitro endpoint is never a QT interval.
"""
from __future__ import annotations

from dataclasses import dataclass, field, replace
from typing import Any, Literal, Sequence

import numpy as np
import pandas as pd

Domain = Literal["clinical", "invitro"]


class ProvenanceError(RuntimeError):
    """Raised when an operation would discard or misuse correction provenance."""


class DomainError(RuntimeError):
    """Raised when an in-vitro endpoint is used as if it were a clinical QT interval."""


@dataclass(frozen=True)
class CorrectionRecord:
    """One applied correction step."""

    step: str                       # e.g. "rate_correction", "baseline_subtraction"
    method: str                     # e.g. "fridericia", "time_matched_predose"
    formula: str                    # human-readable formula actually applied
    parameters: dict[str, Any] = field(default_factory=dict)
    citation: str | None = None
    note: str | None = None

    def __str__(self) -> str:  # pragma: no cover - cosmetic
        p = ", ".join(f"{k}={v!r}" for k, v in self.parameters.items())
        return f"{self.step}:{self.method}({p})"


@dataclass(frozen=True)
class Endpoint:
    """A derived endpoint plus the full cascade that produced it.

    Parameters
    ----------
    values
        The derived numbers.
    name
        Short endpoint name, e.g. ``"QTcF"``, ``"ddQTcF"``, ``"FPDc"``, ``"ddFPDc"``.
    units
        Units of ``values``. ``"ms"`` for intervals, ``"fraction"``, ``"percent"``, ...
    domain
        ``"clinical"`` for patient ECG endpoints, ``"invitro"`` for dish endpoints.
    provenance
        Ordered list of the corrections applied, oldest first.
    """

    values: pd.Series
    name: str
    units: str
    domain: Domain
    provenance: tuple[CorrectionRecord, ...] = ()

    # ---- construction -------------------------------------------------
    def derive(
        self,
        values: pd.Series,
        record: CorrectionRecord,
        *,
        name: str | None = None,
        units: str | None = None,
    ) -> "Endpoint":
        """Return a new Endpoint with ``record`` appended to the cascade."""
        return replace(
            self,
            values=values,
            name=name or self.name,
            units=units or self.units,
            provenance=self.provenance + (record,),
        )

    def subset(self, mask) -> "Endpoint":
        """Select a subset of records, keeping the correction cascade intact.

        ``mask`` is a boolean array or Series aligned with ``values``. Use this
        rather than slicing ``endpoint.values`` directly, which would strip the
        provenance.
        """
        m = np.asarray(mask, dtype=bool)
        if m.size != len(self.values):
            raise ProvenanceError(
                f"mask has {m.size} elements but the endpoint has {len(self.values)}"
            )
        return replace(self, values=self.values[m])

    # ---- introspection ------------------------------------------------
    @property
    def steps(self) -> list[str]:
        return [r.step for r in self.provenance]

    def has_step(self, step: str) -> bool:
        return step in self.steps

    def describe(self) -> str:
        """Multi-line description of the full correction cascade."""
        head = f"{self.name} [{self.units}, {self.domain}] n={len(self.values)}"
        if not self.provenance:
            return head + "\n  (no corrections applied)"
        lines = [head]
        for i, r in enumerate(self.provenance, 1):
            lines.append(f"  {i}. {r.step}: {r.method}")
            lines.append(f"     formula: {r.formula}")
            if r.parameters:
                lines.append(f"     params : {r.parameters}")
            if r.citation:
                lines.append(f"     source : {r.citation}")
            if r.note:
                lines.append(f"     note   : {r.note}")
        return "\n".join(lines)

    def to_frame(self) -> pd.DataFrame:
        """Values as a one-column frame, with the cascade in ``df.attrs``."""
        out = self.values.rename(self.name).to_frame()
        out.attrs["endpoint"] = self.name
        out.attrs["units"] = self.units
        out.attrs["domain"] = self.domain
        out.attrs["provenance"] = [str(r) for r in self.provenance]
        return out

    # ---- guard rails ---------------------------------------------------
    def as_qt_milliseconds(self) -> pd.Series:
        """Return the endpoint as a clinical QT interval in milliseconds.

        Refuses for in-vitro endpoints. Field potential duration and action
        potential duration recorded from a cardiomyocyte monolayer are not QT
        intervals: a monolayer has no ventricular mass, no transmural gradient and
        no conduction geometry, and hiPSC-CMs are electrophysiologically immature.
        Compare the two on a shared relative scale instead - see
        :func:`dish2ecg.bridge.to_shared_scale`.
        """
        if self.domain == "invitro":
            raise DomainError(
                f"{self.name!r} is an in-vitro endpoint and cannot be expressed as a QT "
                "interval in milliseconds. Use percent change or a standardized effect "
                "and compare concentration-response slopes "
                "(dish2ecg.bridge.to_shared_scale)."
            )
        if self.units != "ms":
            raise ProvenanceError(f"{self.name!r} is in {self.units!r}, not milliseconds.")
        return self.values

    def require_steps(self, *steps: str) -> "Endpoint":
        """Assert that the named correction steps were applied; return self."""
        missing = [s for s in steps if s not in self.steps]
        if missing:
            raise ProvenanceError(
                f"{self.name!r} is missing required correction step(s): {missing}. "
                f"Applied: {self.steps}"
            )
        return self


def provenance_columns(endpoint: Endpoint) -> dict[str, str]:
    """Flat dict suitable for writing alongside a results table."""
    return {
        "endpoint": endpoint.name,
        "units": endpoint.units,
        "domain": endpoint.domain,
        "correction_cascade": " | ".join(str(r) for r in endpoint.provenance),
        "citations": "; ".join(
            sorted({r.citation for r in endpoint.provenance if r.citation})
        ),
    }
