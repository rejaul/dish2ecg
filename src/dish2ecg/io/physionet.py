"""Readers for the two public clinical crossover ECG corpora on PhysioNet.

ECGRDVQ  (SCR-002) - Johannesen L et al., Clin Pharmacol Ther 2014;96:549-558.
                     doi:10.1038/clpt.2014.155
                     Dofetilide, quinidine, ranolazine, verapamil, placebo.
ECGDMMLD (SCR-003) - Johannesen L et al., Clin Pharmacol Ther 2015;99:214-223.
                     doi:10.1002/cpt.205
                     Dofetilide alone and with mexiletine or lidocaine;
                     moxifloxacin + diltiazem; placebo.

The package never redistributes these data. ``load_physionet_clinical`` downloads
them on demand and caches them locally.
"""
from __future__ import annotations

import hashlib
import os
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

__all__ = ["PHYSIONET_STUDIES", "PhysionetStudy", "load_physionet_clinical",
           "load_ecgcipa", "ecgcipa_wide", "default_cache_dir"]

_BASE = "https://physionet.org/files"


@dataclass(frozen=True)
class PhysionetStudy:
    slug: str
    version: str
    clinical_csv: str
    description_txt: str
    arm_col: str
    treatments: tuple[str, ...]
    placebo_arm: str
    concentration_layout: str
    citation: str

    @property
    def base_url(self) -> str:
        return f"{_BASE}/{self.slug}/{self.version}"


PHYSIONET_STUDIES: dict[str, PhysionetStudy] = {
    "ecgrdvq": PhysionetStudy(
        slug="ecgrdvq", version="1.0.0",
        clinical_csv="SCR-002.Clinical.Data.csv",
        description_txt="SCR-002.Clinical.Data.Description.txt",
        arm_col="EXTRT",
        treatments=("Dofetilide", "Quinidine Sulph", "Ranolazine",
                    "Verapamil HCL", "Placebo"),
        placebo_arm="Placebo",
        concentration_layout="long",   # PCTEST / PCSTRESN / PCSTRESU
        citation=("Johannesen L et al., Clin Pharmacol Ther 2014;96:549-558. "
                  "doi:10.1038/clpt.2014.155; PhysioNet ECGRDVQ v1.0.0"),
    ),
    "ecgdmmld": PhysionetStudy(
        slug="ecgdmmld", version="1.0.0",
        clinical_csv="SCR-003.Clinical.Data.csv",
        description_txt="SCR-003.Clinical.Data.Description.txt",
        arm_col="TRTA",
        treatments=("Dofetilide", "Lidocaine + Dofetilide", "Mexiletine + Dofetilide",
                    "Moxifloxacin + Diltiazem", "Placebo"),
        placebo_arm="Placebo",
        concentration_layout="wide",   # DOF / LIDO / MEXI / MOXI / MOXI.M2 / DILT
        citation=("Johannesen L et al., Clin Pharmacol Ther 2015;99:214-223. "
                  "doi:10.1002/cpt.205; PhysioNet ECGDMMLD v1.0.0"),
    ),
    "ecgcipa": PhysionetStudy(
        slug="ecgcipa", version="1.0.0",
        clinical_csv="adeg.csv",
        description_txt="ecgcipa_metadata.docx",
        arm_col="TRTA",
        treatments=("Ranolazine", "Verapamil", "Lopinavir+Ritonavir", "Chloroquine",
                    "Dofetilide", "Diltiazem+Dofetilide", "Placebo"),
        placebo_arm="Placebo",
        concentration_layout="adam",   # separate adpc.csv, one row per analyte
        citation=("Vicente J et al., Clin Pharmacol Ther 2019;105:943-953. "
                  "doi:10.1002/cpt.1303; PhysioNet ECGCIPA v1.0.0"),
    ),
}

#: Interval columns in both studies, in milliseconds.
INTERVAL_COLUMNS = ("RR", "PR", "QT", "QRS", "JTPEAK", "TPEAKTEND")


def default_cache_dir() -> Path:
    env = os.environ.get("DISH2ECG_CACHE")
    if env:
        return Path(env)
    return Path.home() / ".cache" / "dish2ecg"


def _download(url: str, dest: Path, timeout: int = 180) -> Path:
    import requests

    dest.parent.mkdir(parents=True, exist_ok=True)
    r = requests.get(url, timeout=timeout)
    r.raise_for_status()
    dest.write_bytes(r.content)
    return dest


def _sha256(p: Path) -> str:
    h = hashlib.sha256()
    h.update(p.read_bytes())
    return h.hexdigest()


def load_ecgcipa(*, cache_dir: str | Path | None = None, download: bool = True
                 ) -> dict[str, pd.DataFrame]:
    """Load the CiPA ECG validation study, which ships as CDISC ADaM tables.

    Unlike ECGRDVQ and ECGDMMLD, this study publishes four analysis datasets
    rather than one flat table:

    ``adsl``  subject level, one row per subject, with the period-1 and period-2
              treatments (``TRT01A`` / ``TRT02A``)
    ``adeg``  ECG analysis data in LONG form - one row per subject, timepoint,
              replicate and parameter (``PARAMCD`` in HR, RR, PR, QRS, QT, QTCF,
              JTP, JTPC, TPTE), with the study's own derived ``BASE`` and ``CHG``
    ``adpc``  pharmacokinetics, one row per analyte per timepoint
    ``addm``  dosing

    Design note: five arms (ranolazine, verapamil, lopinavir+ritonavir,
    chloroquine, placebo) are **parallel**, 10 subjects each; dofetilide and
    diltiazem+dofetilide form a 2-period crossover in a further 10. Placebo
    subtraction is therefore unpaired for the parallel arms - pass
    ``paired=False`` to :func:`dish2ecg.subtract_control`.
    """
    st = PHYSIONET_STUDIES["ecgcipa"]
    cache = Path(cache_dir) if cache_dir else default_cache_dir()
    out: dict[str, pd.DataFrame] = {}
    for name in ("adsl", "adeg", "adpc", "addm"):
        local = cache / st.slug / st.version / f"{name}.csv"
        if not local.exists():
            if not download:
                raise FileNotFoundError(f"{local} not present and download=False")
            _download(f"{st.base_url}/{name}.csv", local)
        out[name] = pd.read_csv(local, low_memory=False)
        out[name].attrs.update(study="ecgcipa", citation=st.citation,
                               sha256=_sha256(local))
    return out


def ecgcipa_wide(adeg: pd.DataFrame) -> pd.DataFrame:
    """Pivot the long ADEG table to one row per ECG with a column per parameter.

    Two of the natural index keys are sparse and must be filled before pivoting,
    or pandas silently drops every row: ``ABLFL`` is 'Y' on the 630 derived
    baseline rows and missing elsewhere, and ``EGREPNUM`` is missing on exactly
    those derived rows (they are not one of the three replicates). Filling them
    with 'N' and 0 keeps all records and makes the baseline flag usable as a
    boolean.
    """
    d = adeg.copy()
    if "ABLFL" in d.columns:
        d["ABLFL"] = d["ABLFL"].astype("object").where(d["ABLFL"].notna(), "N")
    if "EGREPNUM" in d.columns:
        d["EGREPNUM"] = d["EGREPNUM"].fillna(0).astype(int)
    key = [k for k in ["USUBJID", "TRTA", "APERIOD", "ATPTN", "ATPT", "EGREPNUM",
                       "ABLFL"] if k in d.columns]
    w = d.pivot_table(index=key, columns="PARAMCD", values="AVAL",
                      aggfunc="first").reset_index()
    w.columns.name = None
    if len(w) == 0:
        raise RuntimeError("pivot produced no rows; check the index keys for NaNs")
    return w


def load_physionet_clinical(
    study: str, *, cache_dir: str | Path | None = None, download: bool = True,
    verify: bool = True,
) -> pd.DataFrame:
    """Load one clinical table, downloading and caching it if needed.

    The returned frame carries the study metadata in ``df.attrs``: ``study``,
    ``arm_col``, ``placebo_arm``, ``concentration_layout``, ``citation``,
    ``sha256``.
    """
    key = study.lower()
    try:
        st = PHYSIONET_STUDIES[key]
    except KeyError:
        raise ValueError(
            f"unknown study {study!r}; available: {sorted(PHYSIONET_STUDIES)}"
        ) from None

    cache = Path(cache_dir) if cache_dir else default_cache_dir()
    local = cache / st.slug / st.version / st.clinical_csv
    if not local.exists():
        if not download:
            raise FileNotFoundError(
                f"{local} not present and download=False. Fetch it from "
                f"{st.base_url}/{st.clinical_csv}"
            )
        _download(f"{st.base_url}/{st.clinical_csv}", local)

    sums = cache / st.slug / st.version / "SHA256SUMS.txt"
    checksum_status = "not checked"
    if verify:
        if not sums.exists() and download:
            try:
                _download(f"{st.base_url}/SHA256SUMS.txt", sums)
            except Exception:
                pass
        if sums.exists():
            want = None
            for line in sums.read_text(errors="replace").splitlines():
                parts = line.split()
                if len(parts) >= 2 and parts[-1].endswith(st.clinical_csv):
                    want = parts[0]
                    break
            got = _sha256(local)
            checksum_status = ("ok" if want and want.lower() == got.lower()
                               else f"MISMATCH (expected {want}, got {got})"
                               if want else "no entry in SHA256SUMS.txt")
            if want and want.lower() != got.lower():
                raise RuntimeError(
                    f"checksum mismatch for {local}: expected {want}, got {got}"
                )

    df = pd.read_csv(local)
    df.attrs.update(
        study=key, arm_col=st.arm_col, placebo_arm=st.placebo_arm,
        treatments=list(st.treatments),
        concentration_layout=st.concentration_layout,
        citation=st.citation, source_url=f"{st.base_url}/{st.clinical_csv}",
        sha256=_sha256(local), checksum_status=checksum_status,
    )
    return df


def concentration_series(df: pd.DataFrame, analyte: str) -> pd.Series:
    """Plasma concentration of ``analyte`` for every row, NaN where not measured.

    Handles both layouts: ECGRDVQ stores one analyte per row in PCTEST/PCSTRESN,
    ECGDMMLD stores one column per analyte.
    """
    layout = df.attrs.get("concentration_layout")
    if layout == "long":
        if not {"PCTEST", "PCSTRESN"} <= set(df.columns):
            raise KeyError("expected PCTEST and PCSTRESN columns for the long layout")
        return df["PCSTRESN"].where(df["PCTEST"] == analyte).astype(float)
    if layout == "wide":
        if analyte not in df.columns:
            raise KeyError(
                f"{analyte!r} not among the concentration columns; "
                f"candidates: {[c for c in df.columns if c.isupper()]}"
            )
        return pd.to_numeric(df[analyte], errors="coerce")
    raise ValueError(
        "df.attrs['concentration_layout'] is missing; load the table with "
        "load_physionet_clinical()"
    )
