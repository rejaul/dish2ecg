#!/usr/bin/env python
"""Join the published CiPA myocyte workbook to the concentration lookup.

    python scripts/build_wells.py \
        --workbook Blinova_etal_2018_data.xlsx \
        --lookup data/cipa_concentration_lookup.csv \
        --out data/cipa_myocyte_wells.csv

The workbook is the well-level release of the FDA/HESI CiPA myocyte validation
study, distributed from the CiPA Data Resources page. It is **not** redistributed
with this repository. It carries a concentration *index* (1-4) but not the molar
concentration; the lookup supplies the molar concentration, the multiple of free
Cmax, and the Cmax itself, read from the per-drug tables printed under the figure
panels of the open-access supplement (Data S1-S3 of
doi:10.1016/j.celrep.2018.08.079).
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

#: Documented in Supplemental Table 1 of the source paper.
CELL_LINE = {"CDI": "iCell2", "AXG": "Cor.4U"}
PLATFORM = {"AXN": "Maestro (Axion Biosystems)",
            "MCS": "MEA2100 (Multichannel Systems)",
            "CLY": "CellOPTIQ (Clyde Biosciences)",
            "ECR": "CardioECR (ACEA Biosciences)",
            "AMD": "AlphaMED64 (Alpha MED Scientific)",
            "ACA": "CardioECR (ACEA Biosciences)"}
RISK = {"H": "high", "M": "intermediate", "L": "low"}
#: One compound appears under two spellings in the workbook.
DRUG_FIX = {"D,l,Sotalol": "D,l Sotalol"}


def main(args) -> None:
    w = pd.read_excel(args.workbook, sheet_name=0)
    n_fixed = int(w["Drug_Name"].isin(DRUG_FIX).sum())
    w["drug"] = w["Drug_Name"].replace(DRUG_FIX)
    w["drug_name_normalised"] = w["Drug_Name"].isin(DRUG_FIX)
    w["cell_line"] = w["Cell_type"].map(CELL_LINE)
    w["platform_name"] = w["Platform"].map(PLATFORM)
    w["risk_category"] = w["risk"].map(RISK)

    lut = pd.read_csv(args.lookup)
    out = w.merge(lut.rename(columns={"conc_index": "conc"}), on=["drug", "conc"],
                  how="left", validate="many_to_one")

    missing = int(out["conc_molar"].isna().sum())
    print(f"rows: {len(out)}; drug-name spellings normalised: {n_fixed}; "
          f"rows without a concentration lookup: {missing}")
    print("cell lines:", out["cell_line"].value_counts().to_dict())
    print("platforms :", out["platform_name"].value_counts().to_dict())
    if missing:
        print("UNMATCHED drugs:",
              sorted(out.loc[out["conc_molar"].isna(), "drug"].unique()))
    out.to_csv(args.out, index=False)
    print("wrote", args.out)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--workbook", type=Path, required=True)
    ap.add_argument("--lookup", type=Path,
                    default=Path("data/cipa_concentration_lookup.csv"))
    ap.add_argument("--out", type=Path, default=Path("data/cipa_myocyte_wells.csv"))
    main(ap.parse_args())
