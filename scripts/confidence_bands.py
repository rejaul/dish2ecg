"""Concentration-response with confidence bands, one panel per arm.

Each panel is read at the geometric-mean OBSERVED clinical free Cmax and nowhere
else, so that point is marked with its interval in every panel. Bands use
`dish2ecg.viz.plot_confidence_band`, which reproduces the interval the model
layer reports (t quantile on n_groups - 1), so the band and the quoted effect
cannot disagree.

Drug: dofetilide, the one compound with a clear signal in all three arms.

Writes results/confidence_bands.png and results/confidence_bands.csv.

    python scripts/confidence_bands.py
"""
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from dish2ecg import viz

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
PNG = RES / "confidence_bands.png"
CSV = RES / "confidence_bands.csv"

DISH, CLIN, SIL, RED, GREY = "#C8541A", "#2B5D9B", "#4A4A4A", "#B3242B", "#8A8A8A"
DRUG = "Dofetilide"
LEVEL = 0.90


def main():
    mpl.rcParams.update({"font.family": "sans-serif", "font.size": 8, "figure.dpi": 200,
                         "figure.facecolor": "white", "axes.facecolor": "white",
                         "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(1, 3, figsize=(10.6, 3.5))
    rows = []

    # ---------------- in vitro
    iv = pd.read_csv(RES / "invitro_concentration_response.csv")
    r = iv[(iv.drug == DRUG) & (iv.analysis == "all_concentrations") & (iv.endpoint == "ddFPDc")].iloc[0]
    bridge = pd.read_csv(RES / "dish_vs_clinic_bridge.csv")
    match = float(bridge.loc[bridge.drug == DRUG, "match_point_folds_free_cmax"].iloc[0])
    ax = axes[0]
    xs = np.linspace(-1.0, 0.8, 200)
    viz.plot_confidence_band(xs, slope=r.slope, intercept=r.intercept, slope_se=r.slope_se,
                             intercept_se=r.intercept_se, centre=float(r.at_concentration),
                             df=int(r.n_groups) - 1, level=LEVEL, ax=ax, color=DISH,
                             mark_centre=False)
    f, lo, hi = viz.confidence_band([np.log10(match)], slope=r.slope, intercept=r.intercept,
                                    slope_se=r.slope_se, intercept_se=r.intercept_se,
                                    centre=float(r.at_concentration), df=int(r.n_groups) - 1,
                                    level=LEVEL)
    ax.axvline(np.log10(match), color=GREY, lw=0.7, ls=(0, (2, 2)), zorder=0)
    ax.errorbar([np.log10(match)], f, yerr=[[f[0] - lo[0]], [hi[0] - f[0]]], fmt="o", ms=4.5,
                color=DISH, capsize=2.4, lw=1.1, zorder=5, mfc="white", mew=1.1)
    ax.annotate(f"{f[0]:+.0f} ms\n[{lo[0]:.0f}, {hi[0]:.0f}]", xy=(np.log10(match), f[0]),
                xytext=(8, -4), textcoords="offset points", fontsize=6.6, color=DISH, va="top")
    ax.axhline(0, color="#CCCCCC", lw=0.6)
    ax.set_xlabel("log$_{10}$(\u00d7 free C$_{max}$)")
    ax.set_ylabel("\u0394\u0394FPDc (ms)")
    ax.set_title(f"a   In vitro \u00b7 {int(r.n_obs)} wells, {int(r.n_groups)} datasets", loc="left", fontsize=8.4)
    ax.text(0.03, 0.03, "no in-vitro equivalent\nof the 10 ms rule", transform=ax.transAxes,
            fontsize=6.0, color=GREY, va="bottom", linespacing=1.3)
    rows.append(dict(arm="in vitro", endpoint="ddFPDc", drug=DRUG, at=f"{match:.3f}x free Cmax",
                     effect=f[0], lo=lo[0], hi=hi[0], df=int(r.n_groups) - 1))

    # ---------------- in vivo
    cq = pd.read_csv(RES / "concentration_qtc_results.csv")
    q = cq[(cq.endpoint == "ddQTcF") & (cq.study == "ecgrdvq") & (cq.arm == DRUG)].iloc[0]
    clin = pd.read_csv(RES / "clinical_ddqtc.csv")
    pts = clin[(clin.study == "ecgrdvq") & (clin.arm == DRUG) & (~clin.is_baseline)]
    ax = axes[1]
    ax.scatter(pts.conc_ngml, pts.ddQTcF, s=3, color=CLIN, alpha=0.25, lw=0)
    xs = np.linspace(0, float(pts.conc_ngml.max()), 200)
    viz.plot_confidence_band(xs, slope=q.slope, intercept=q.intercept, slope_se=q.slope_se,
                             intercept_se=q.intercept_se, centre=float(q.at_concentration),
                             df=int(q.n_groups) - 1, level=LEVEL, ax=ax, color=CLIN,
                             centre_label=f"{q.effect_at_concentration:+.0f} ms\n"
                                          f"[{q.effect_ci90_lower:.0f}, {q.effect_ci90_upper:.0f}]")
    ax.axhline(0, color="#CCCCCC", lw=0.6)
    ax.axhline(10, color=RED, lw=0.9, ls=(0, (4, 2)))
    ax.text(0.98, 10, "10 ms, ICH E14 ", transform=ax.get_yaxis_transform(), ha="right",
            va="bottom", fontsize=6.2, color=RED)
    ax.set_xlabel("plasma concentration (ng/mL)")
    ax.set_ylabel("\u0394\u0394QTcF (ms)")
    ax.set_title(f"b   In vivo \u00b7 {int(q.n_obs)} records, {int(q.n_groups)} subjects", loc="left", fontsize=8.4)
    rows.append(dict(arm="in vivo", endpoint="ddQTcF", drug=DRUG,
                     at=f"{q.at_concentration:.2f} ng/mL (geometric-mean observed Cmax)",
                     effect=q.effect_at_concentration, lo=q.effect_ci90_lower,
                     hi=q.effect_ci90_upper, df=int(q.n_groups) - 1))

    # ---------------- in silico
    od = pd.read_csv(RES / "ord_apd90.csv")
    o = od[(od.drug == DRUG.lower()) & od.apd90.notna() & (od.v_rest < -60)]
    wide = o.pivot_table(index="multiple", columns="cell", values="dAPD90").dropna()
    ax = axes[2]
    ax.fill_between(wide.index, wide.min(axis=1), wide.max(axis=1), color=SIL, alpha=0.17, lw=0)
    ax.plot(wide.index, wide["mid"], color=SIL, lw=1.4, marker="o", ms=2.4, mfc=SIL, mec="none")
    at = float(np.interp(match, wide.index, wide["mid"]))
    band = (float(np.interp(match, wide.index, wide.min(axis=1))),
            float(np.interp(match, wide.index, wide.max(axis=1))))
    ax.axvline(match, color=GREY, lw=0.7, ls=(0, (2, 2)), zorder=0)
    ax.errorbar([match], [at], yerr=[[at - band[0]], [band[1] - at]], fmt="o", ms=4.5, color=SIL,
                capsize=2.4, lw=1.1, zorder=5, mfc="white", mew=1.1)
    ax.annotate(f"{at:+.0f} ms\n[{band[0]:.0f}, {band[1]:.0f}]", xy=(match, at), xytext=(8, -4),
                textcoords="offset points", fontsize=6.6, color=SIL, va="top")
    ax.axhline(0, color="#CCCCCC", lw=0.6)
    ax.set_xscale("log")
    ax.set_xlabel("\u00d7 free C$_{max}$")
    ax.set_ylabel("\u0394APD$_{90}$ (ms)")
    ax.set_title("c   In silico \u00b7 mid-myocardial cell", loc="left", fontsize=8.4)
    ax.text(0.03, 0.95, "band = range over the three cell\ntypes, NOT a confidence interval:\n"
                        "the simulation has no sampling error", transform=ax.transAxes,
            fontsize=6.0, color=GREY, va="top", linespacing=1.35)
    rows.append(dict(arm="in silico", endpoint="dAPD90", drug=DRUG, at=f"{match:.3f}x free Cmax",
                     effect=at, lo=band[0], hi=band[1], df=np.nan))

    out = pd.DataFrame(rows)
    out["interval_meaning"] = ["90 % CI, t on n_groups-1", "90 % CI, t on n_groups-1",
                               "range over endo/epi/mid, not a CI"]
    out.to_csv(CSV, index=False)
    fig.subplots_adjust(left=0.065, right=0.99, bottom=0.165, top=0.88, wspace=0.33)
    fig.savefig(PNG, dpi=200)
    print(out.to_string(index=False))
    print(f"\nwrote {PNG}\nwrote {CSV}")


if __name__ == "__main__":
    main()
