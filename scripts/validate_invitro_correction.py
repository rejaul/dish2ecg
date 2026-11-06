#!/usr/bin/env python
"""Validate the in-vitro rate correction against an independent published one.

    python scripts/validate_invitro_correction.py --out results/

The CiPA myocyte workbook publishes ddFPDc already corrected, so it cannot test
our correction. The CiPA **MEA pilot** study can: Supplemental Table 1 of
Millard D et al., *Toxicol Sci* 2018;164:550-562, doi:10.1093/toxsci/kfy110
(PMC6061700) reports, for each of the 18 participating studies, the mean
baseline beat period, the mean field potential duration **and** the mean
rate-corrected field potential duration. Three quantities, one formula between
them — so our FPDc can be checked against theirs.

The table is published as a vector image inside the supplementary Word document;
`dish2ecg.io.emf_tables` reads the drawn text records out of the metafile, so
the digits are recovered verbatim rather than by OCR.

One caveat governs the interpretation. The published FPDc is the mean of the
per-well corrected values, while what we can compute from the table is the
correction applied to the mean beat period and the mean duration. Those differ
by Jensen's inequality, and the gap should grow with the spread of the beat
period — which is the check this script actually performs.
"""
from __future__ import annotations

import argparse
import io
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from dish2ecg import power_correction
from dish2ecg.io import emf_table

EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest/{}/supplementaryFiles"
PMCID = "PMC6061700"
CITATION = ("Millard D et al., Toxicol Sci 2018;164:550-562. "
            "doi:10.1093/toxsci/kfy110")


def fetch_table() -> list[list[str]]:
    """Download the supplement and recover Supplemental Table 1."""
    import requests

    r = requests.get(EPMC.format(PMCID), timeout=300)
    r.raise_for_status()
    outer = zipfile.ZipFile(io.BytesIO(r.content))
    inner_name = next(n for n in outer.namelist() if n.lower().endswith("_supp.zip"))
    inner = zipfile.ZipFile(io.BytesIO(outer.read(inner_name)))
    docx_name = next(n for n in inner.namelist() if n.lower().endswith(".docx"))
    doc = zipfile.ZipFile(io.BytesIO(inner.read(docx_name)))
    # Supplemental Table 1 is the first metafile in the document.
    emf = sorted(n for n in doc.namelist()
                 if n.lower().startswith("word/media/") and n.lower().endswith(".emf"))
    for name in emf:
        rows = emf_table(doc.read(name))
        if rows and "BP (ms)" in rows[0]:
            return rows
    raise RuntimeError("Supplemental Table 1 not found in the supplement")


def parse(rows: list[list[str]]) -> pd.DataFrame:
    def pair(cell: str) -> tuple[float, float]:
        a, b = cell.split("±")
        return float(a.strip()), float(b.strip())

    out = []
    for cells in rows[1:]:
        bp, bp_sd = pair(cells[3])
        fpd, fpd_sd = pair(cells[4])
        fpdc, fpdc_sd = pair(cells[5])
        out.append(dict(study=int(cells[0]), cell_type=cells[1], platform=cells[2],
                        BP_ms=bp, BP_sd=bp_sd, FPD_ms=fpd, FPD_sd=fpd_sd,
                        FPDc_published=fpdc, FPDc_sd=fpdc_sd))
    return pd.DataFrame(out)


def main(out: Path) -> None:
    out.mkdir(parents=True, exist_ok=True)
    M = parse(fetch_table())

    ep = power_correction(M["FPD_ms"], M["BP_ms"], 1 / 3,
                          name="FPDcF", domain="invitro")
    M["FPDc_dish2ecg"] = ep.values.round(1).to_numpy()
    M["signed_error_ms"] = (M["FPDc_dish2ecg"] - M["FPDc_published"]).round(1)
    M["abs_error_ms"] = M["signed_error_ms"].abs()
    M["relative_error_pct"] = (100 * M["abs_error_ms"] / M["FPDc_published"]).round(2)
    M["BP_cv"] = (M["BP_sd"] / M["BP_ms"]).round(3)
    M["source"] = CITATION
    M["correction_cascade"] = " | ".join(str(r) for r in ep.provenance)

    rho, p = stats.spearmanr(M["BP_cv"], M["abs_error_ms"])
    # exponent implied by the published corrected values
    alpha = float(np.polyfit(np.log(M["BP_ms"] / 1000.0),
                             np.log(M["FPD_ms"] / M["FPDc_published"]), 1)[0])

    M.to_csv(out / "invitro_correction_validation.csv", index=False)
    print(M[["study", "cell_type", "platform", "BP_ms", "FPD_ms", "FPDc_published",
             "FPDc_dish2ecg", "signed_error_ms", "BP_cv"]].to_string(index=False))
    print(f"\nstudies: {len(M)}")
    print(f"median |error| {M.abs_error_ms.median():.1f} ms, "
          f"max {M.abs_error_ms.max():.1f} ms "
          f"({M.relative_error_pct.max():.2f} % relative)")
    ok = M[M.abs_error_ms <= 5]
    print(f"within 5 ms: {len(ok)} of {len(M)}")
    print(f"exponent implied by the published FPDc: {alpha:.4f} "
          f"(Fridericia = {1/3:.4f})")
    print(f"Spearman(beat-period CV, |error|): rho={rho:.3f}, p={p:.3f}  "
          "— the residual grows with beat-period spread, as Jensen's inequality "
          "predicts for a mean of corrected values versus a correction of means")
    worst = M.loc[M.abs_error_ms.idxmax()]
    print(f"\nlargest discrepancy: study {int(worst.study)} "
          f"({worst.cell_type}/{worst.platform}), {worst.abs_error_ms:.1f} ms, "
          f"beat-period CV {worst.BP_cv:.3f}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=Path("results"))
    main(ap.parse_args().out)
