"""Does a relative scale (percent change, or ratio to a reference drug) make the dish and the
clinic agree, where milliseconds do not?

Motivation: FPD and QT are different observables, so agreement in milliseconds is not defined.
The natural repair is to normalise - express each effect as a percent of its own baseline, or as
a ratio to a reference drug measured in the same system. This script tests whether that works.

Writes results/relative_scale_test.csv.

    python scripts/relative_scale_test.py
"""
import re
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
OUT = RES / "relative_scale_test.csv"

# ECGRDVQ arm labels for the four bridge drugs
ARM = {"Dofetilide": "Dofetilide", "Quinidine": "Quinidine Sulph",
       "Ranolazine": "Ranolazine", "Verapamil": "Verapamil HCL"}
REFERENCE = "Dofetilide"   # the reference drug for the ratio normalisation


def _ci(s):
    return [float(x) for x in re.findall(r"-?\d+\.?\d*", str(s))]


def main():
    bridge = pd.read_csv(RES / "dish_vs_clinic_bridge.csv")
    clin = pd.read_csv(RES / "clinical_ddqtc.csv")
    e = clin[clin.study == "ecgrdvq"]
    # clinical baseline QTcF: the drug-free pre-dose level in that arm
    baseline_qtcf = {d: float(e[(e.arm == a) & (e.is_baseline)].QTcF.mean()) for d, a in ARM.items()}

    # The CiPA myocyte workbook ships only ddFPDc, with no baseline level, so a percent change
    # cannot be formed from it. We import an EXTERNAL baseline: the median control FPDc over the
    # 18 Millard 2018 studies, which use the same two cell lines and the same platforms.
    millard = pd.read_csv(RES / "invitro_correction_validation.csv")
    baseline_fpdc = float(millard.FPDc_published.median())

    ref_clin = float(bridge.loc[bridge.drug == REFERENCE, "clin_ddQTcF_ms"].iloc[0])
    ref_dish = float(bridge.loc[bridge.drug == REFERENCE, "dish_ddFPDc_ms"].iloc[0])

    rows = []
    for _, r in bridge.iterrows():
        lo, hi = _ci(r.clin_ci)[:2]
        clin_pct = 100 * r.clin_ddQTcF_ms / baseline_qtcf[r.drug]
        dish_pct = 100 * r.dish_ddFPDc_ms / baseline_fpdc
        rows.append(dict(
            drug=r.drug, risk=r.risk,
            clin_ddQTcF_ms=r.clin_ddQTcF_ms, dish_ddFPDc_ms=r.dish_ddFPDc_ms,
            clin_baseline_ms=round(baseline_qtcf[r.drug], 1), dish_baseline_ms=round(baseline_fpdc, 1),
            clin_pct_of_baseline=round(clin_pct, 3), dish_pct_of_baseline=round(dish_pct, 3),
            ratio_raw_ms=round(r.dish_ddFPDc_ms / r.clin_ddQTcF_ms, 3),
            ratio_percent_scale=round(dish_pct / clin_pct, 3),
            clin_rel_to_reference=round(r.clin_ddQTcF_ms / ref_clin, 3),
            dish_rel_to_reference=round(r.dish_ddFPDc_ms / ref_dish, 3),
            ratio_reference_scale=round((r.dish_ddFPDc_ms / ref_dish) / (r.clin_ddQTcF_ms / ref_clin), 3),
            clin_ci_spans_zero=bool(lo < 0 < hi),
        ))
    out = pd.DataFrame(rows)
    out["reference_drug"] = REFERENCE
    out["dish_baseline_source"] = "median control FPDc over 18 Millard 2018 studies (external)"
    out.to_csv(OUT, index=False)

    # assertions that keep the conclusions honest
    pure = out[out.drug.isin(["Dofetilide", "Quinidine"])]
    assert abs(pure.ratio_raw_ms.iloc[0] - pure.ratio_raw_ms.iloc[1]) < 0.01, \
        "the two pure hERG blockers no longer share one dish/clinic scale factor"
    vera = out[out.drug == "Verapamil"].iloc[0]
    assert vera.clin_ci_spans_zero, "verapamil's clinical interval no longer spans zero"
    print(out[["drug", "clin_pct_of_baseline", "dish_pct_of_baseline",
               "ratio_raw_ms", "ratio_percent_scale", "ratio_reference_scale"]].to_string(index=False))
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
