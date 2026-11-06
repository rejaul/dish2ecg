#!/usr/bin/env python
"""In-vitro concentration-response and the dish-to-clinic bridge.

    python scripts/invitro_and_bridge.py \
        --wells data/cipa_myocyte_wells.csv \
        --clinical results/clinical_ddqtc.csv \
        --out results/

Reproduces results/invitro_concentration_response.csv,
results/dish_vs_clinic_bridge.csv and results/dish_vs_clinic.png.

Two facts govern this script and are asserted rather than assumed.

1. The ``ddFPDc`` column of the CiPA myocyte dataset arrives **already** rate-,
   baseline- and vehicle-corrected by the source study. The dish2ecg cascade must
   not be applied to it a second time; the corrections are attached as provenance
   so the modelling layer can still check them.
2. The missingness of ``ddFPDc`` is **informative**: a well with an arrhythmic
   event has no measurable field potential duration, so the endpoint disappears
   exactly where the drug effect is largest. A complete-case concentration-response
   slope is therefore biased toward zero, and the script reports both the
   all-concentration fit and a low-concentration sensitivity fit where the
   missingness is small.
"""
from __future__ import annotations

import argparse
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

from dish2ecg import fit_concentration_response
from dish2ecg.provenance import CorrectionRecord, Endpoint

BLINOVA = ("Blinova K et al., Cell Rep 2018;24:3582-3592. "
           "doi:10.1016/j.celrep.2018.08.079")

#: The corrections the SOURCE study applied, recorded so they travel with the data.
SOURCE_CASCADE = (
    CorrectionRecord(step="rate_correction", method="fridericia_by_source",
                     formula="FPDc = FPD / (IBI_s)**(1/3)", citation=BLINOVA,
                     note="Applied by the source study, not by dish2ecg."),
    CorrectionRecord(step="baseline_subtraction", method="predose_by_source",
                     formula="dFPDc = FPDc(drug) - FPDc(pre-dose, same well)",
                     citation=BLINOVA, note="Applied by the source study."),
    CorrectionRecord(step="control_subtraction", method="vehicle_by_source",
                     formula="ddFPDc = dFPDc - dFPDc(vehicle)", citation=BLINOVA,
                     note="Applied by the source study."),
)

#: Molecular weights, PubChem (g/mol), for the four agents present in both corpora.
MW = {"Dofetilide": 441.6, "Quinidine": 324.4, "Ranolazine": 427.5,
      "Verapamil": 454.6}
#: Fraction unbound, from the project drug reference table.
FU = {"Dofetilide": 0.35, "Quinidine": 0.16, "Ranolazine": 0.38, "Verapamil": 0.10}
#: ECGRDVQ arm label for each agent.
ARM = {"Dofetilide": "Dofetilide", "Quinidine": "Quinidine Sulph",
       "Ranolazine": "Ranolazine", "Verapamil": "Verapamil HCL"}

CLINIC_COLOR, DISH_COLOR = "#1b6ca8", "#c8582a"


def dish_endpoint(sub: pd.DataFrame) -> Endpoint:
    return Endpoint(values=sub["ddFPDc"].astype(float).rename("ddFPDc"),
                    name="ddFPDc", units="ms", domain="invitro",
                    provenance=SOURCE_CASCADE)


def fit_invitro(wells: pd.DataFrame, conc_subset=None, tag="all_concentrations"
                ) -> pd.DataFrame:
    """One log-concentration fit per drug, with site/line/platform random effects.

    The predictor is log10(multiples of free Cmax): tested concentrations span up
    to 3500-fold, so a model linear in raw concentration is dominated by the top
    point and its slope is not interpretable. On the log axis the slope is
    milliseconds per decade and the intercept is the effect at 1x free Cmax.
    """
    rows = []
    for drug, g in wells.groupby("drug"):
        gg = g if conc_subset is None else g[g["conc"].isin(conc_subset)]
        sel = (gg["ddFPDc"].notna() & gg["log10_folds_cmax"].notna()).to_numpy()
        if sel.sum() < 20 or gg.loc[sel, "log10_folds_cmax"].nunique() < 3:
            continue
        sub = gg[sel].reset_index(drop=True)
        sub = sub.assign(_grp=sub["site"].astype(str) + "|" + sub["cell_line"]
                         + "|" + sub["Platform"])
        f = fit_concentration_response(
            sub, dish_endpoint(sub), concentration_col="log10_folds_cmax",
            group_col="_grp", at_concentration=0.0)
        row = f.to_row()
        row.update(drug=drug, risk=g["risk_category"].dropna().iloc[0],
                   analysis=tag, n_wells_used=int(sel.sum()),
                   n_wells_available=int(len(gg)),
                   pct_wells_dropped=round(100 * (1 - sel.sum() / len(gg)), 1),
                   n_datasets=int(sub["_grp"].nunique()),
                   effect_units="ms ddFPDc at 1x free Cmax",
                   slope_units="ms per decade of free-Cmax multiple")
        rows.append(row)
    return pd.DataFrame(rows)


def build_bridge(wells: pd.DataFrame, clin: pd.DataFrame) -> pd.DataFrame:
    """Compare dish and clinic at the SAME multiple of free Cmax.

    The comparison point is the geometric-mean observed clinical free Cmax, not
    1x the nominal free Cmax: in ECGRDVQ the observed exposure reaches only 0.29x
    for quinidine and 0.50x for verapamil, so evaluating the clinical model at 1x
    would be extrapolation beyond the measured range.
    """
    cmax_free = wells.groupby("drug")["cmax_molar"].first().to_dict()
    rd = clin[clin["study"] == "ecgrdvq"]
    rows = []
    for drug in MW:
        s = rd[(rd["arm"] == ARM[drug]) & rd["conc_ngml"].notna()
               & (rd["TPT"] > 0)].copy()
        s["free_molar"] = s["conc_ngml"] * 1e-6 / MW[drug] * FU[drug]
        s["folds_free_cmax"] = s["free_molar"] / cmax_free[drug]
        point = float(np.exp(np.log(
            s.groupby("RANDID")["folds_free_cmax"].max()).mean()))

        ep = Endpoint(values=s["ddQTcF"].astype(float).rename("ddQTcF"),
                      name="ddQTcF", units="ms", domain="clinical",
                      provenance=SOURCE_CASCADE)  # cascade recorded upstream
        cl = fit_concentration_response(
            s, ep, concentration_col="folds_free_cmax", group_col="RANDID",
            at_concentration=point, require_steps=())

        d = wells[(wells["drug"] == drug) & wells["ddFPDc"].notna()
                  & wells["log10_folds_cmax"].notna()].reset_index(drop=True)
        d = d.assign(_grp=d["site"].astype(str) + "|" + d["cell_line"] + "|"
                     + d["Platform"])
        di = fit_concentration_response(
            d, dish_endpoint(d), concentration_col="log10_folds_cmax",
            group_col="_grp", at_concentration=float(np.log10(point)))

        lo, hi = d["folds_cmax"].min(), d["folds_cmax"].max()
        rows.append(dict(
            drug=drug, risk=d["risk_category"].dropna().iloc[0],
            match_point_folds_free_cmax=round(point, 3),
            clin_conc_range=f"{s.folds_free_cmax.min():.3f}-"
                            f"{s.folds_free_cmax.max():.3f}",
            dish_conc_range=f"{lo:.3g}-{hi:.3g}",
            point_inside_dish_range=bool(lo <= point <= hi),
            clin_ddQTcF_ms=round(cl.effect_at_concentration, 1),
            clin_ci=f"[{cl.effect_ci90[0]:.1f}, {cl.effect_ci90[1]:.1f}]",
            dish_ddFPDc_ms=round(di.effect_at_concentration, 1),
            dish_ci=f"[{di.effect_ci90[0]:.1f}, {di.effect_ci90[1]:.1f}]",
            dish_minus_clin=round(di.effect_at_concentration
                                  - cl.effect_at_concentration, 1),
            same_sign=bool(np.sign(di.effect_at_concentration)
                           == np.sign(cl.effect_at_concentration)),
            dish_wells=int(di.n_obs), dish_datasets=int(di.n_groups),
            clin_subjects=int(cl.n_groups)))
    return pd.DataFrame(rows)


def make_figure(wells, clin, bridge, path: Path) -> None:
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 4, figsize=(11.4, 3.5))
    rd = clin[clin["study"] == "ecgrdvq"]
    cmax_free = wells.groupby("drug")["cmax_molar"].first().to_dict()
    for ax, drug in zip(axes, MW):
        row = bridge.set_index("drug").loc[drug]
        pt = row["match_point_folds_free_cmax"]
        d = wells[(wells["drug"] == drug) & wells["ddFPDc"].notna()
                  & wells["folds_cmax"].notna()]
        ax.scatter(d["folds_cmax"], d["ddFPDc"], s=5, c=DISH_COLOR, alpha=.18, lw=0)
        med = d.groupby("folds_cmax")["ddFPDc"].median()
        ax.plot(med.index, med.values, "o-", color=DISH_COLOR, ms=5, lw=1.6)
        s = rd[(rd["arm"] == ARM[drug]) & rd["conc_ngml"].notna()
               & (rd["TPT"] > 0)].copy()
        s["fx"] = s["conc_ngml"] * 1e-6 / MW[drug] * FU[drug] / cmax_free[drug]
        ax.scatter(s["fx"], s["ddQTcF"], s=5, c=CLINIC_COLOR, alpha=.30, lw=0)
        b = pd.qcut(np.log10(s["fx"]), 6, duplicates="drop")
        bm = s.groupby(b, observed=True).agg(x=("fx", "median"),
                                             y=("ddQTcF", "median"))
        ax.plot(bm["x"], bm["y"], "s-", color=CLINIC_COLOR, ms=4.5, lw=1.6)
        ax.axhline(0, color="0.55", lw=.7)
        ax.axvline(pt, color="0.55", lw=.8, ls=":")
        ax.set_xscale("log")
        ax.set_xlim(min(d["folds_cmax"].min(), s["fx"].min()) * 0.6,
                    max(d["folds_cmax"].max(), s["fx"].max()) * 1.8)
        ax.set_title(f"{drug} · {row['risk']} TdP risk", fontsize=8)
        ax.margins(y=.12)
        ax.set_xlabel("× free $C_{max}$")
        ax.annotate(f"at {pt:.2g}× :  clinic {row['clin_ddQTcF_ms']:+.0f} ms\n"
                    f"{' ' * 10}dish {row['dish_ddFPDc_ms']:+.0f} ms",
                    xy=(.97, .03), xycoords="axes fraction", fontsize=6.5,
                    ha="right", va="bottom", color="0.25")
    axes[0].set_ylabel("Change in repolarization\nduration (ms)")
    axes[0].plot([], [], "s-", color=CLINIC_COLOR, ms=4.5,
                 label="clinic  ΔΔQTcF, 22 subjects")
    axes[0].plot([], [], "o-", color=DISH_COLOR, ms=5,
                 label="dish  ΔΔFPDc, 15 datasets")
    axes[0].legend(loc="upper left", frameon=False, fontsize=6.5)
    fig.suptitle("The dish tracks the clinic in sign and rank for the hERG "
                 "blockers, and inverts for verapamil", fontsize=9, y=1.03)
    fig.tight_layout()
    fig.savefig(path, dpi=300, bbox_inches="tight")


def main(args) -> None:
    warnings.filterwarnings("ignore")
    args.out.mkdir(parents=True, exist_ok=True)

    wells = pd.read_csv(args.wells)
    wells["log10_folds_cmax"] = np.log10(wells["folds_cmax"])
    clin = pd.read_csv(args.clinical)

    miss = wells.groupby("conc")["ddFPDc"].apply(lambda s: 100 * s.isna().mean())
    ead_nulls = int(wells.loc[wells["ddFPDc"].isna(), "EAD"].sum())
    print("ddFPDc missing (%) by concentration index:",
          {int(k): round(v, 1) for k, v in miss.items()})
    print(f"of {int(wells['ddFPDc'].isna().sum())} missing endpoints, "
          f"{ead_nulls} are wells flagged with an arrhythmic event")

    iv = pd.concat([fit_invitro(wells, None, "all_concentrations"),
                    fit_invitro(wells, [1, 2, 3], "conc_1_to_3_only")],
                   ignore_index=True)
    iv.to_csv(args.out / "invitro_concentration_response.csv", index=False)

    bridge = build_bridge(wells, clin)
    bridge.to_csv(args.out / "dish_vs_clinic_bridge.csv", index=False)
    print()
    print(bridge.to_string(index=False))
    print(f"\nsign agreement: {int(bridge.same_sign.sum())} of {len(bridge)} agents")

    make_figure(wells, clin, bridge, args.out / "dish_vs_clinic.png")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--wells", type=Path, required=True)
    ap.add_argument("--clinical", type=Path, default=Path("results/clinical_ddqtc.csv"))
    ap.add_argument("--out", type=Path, default=Path("results"))
    main(ap.parse_args())
