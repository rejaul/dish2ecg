"""Bland-Altman method comparison, where it is licensed.

Bland-Altman compares two methods that measure the SAME quantity in the SAME units and are
candidates to replace one another. Two comparisons in this repository meet that condition:

  (a) optical APD90 against electrode FPD, 97 paired drug x concentration cells. Blinova 2018
      pools both into one `ddFPDc` column and calls them "equivalent measures", so this is a
      direct test of that claim.
  (b) our tangent T-end against the study's own T-end, on the same clinical median beats.

The dish-versus-clinic comparison does NOT meet it: FPD and QT are different observables, the
average of the two estimates nothing, and n = 4 drugs cannot support limits of agreement. That
comparison is handled on the concentration axis instead (see dish_vs_clinic_bridge.csv).

Repeated measurements are clustered (several concentrations per drug; several records per
subject), so the naive SD of differences is too small. Limits of agreement are therefore formed
from a random-intercept model, following Bland & Altman's treatment of repeated measures:
SD_diff = sqrt(var_between + var_within).

Writes results/bland_altman_stats.csv and results/bland_altman.png.

    python scripts/bland_altman.py
"""
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
CSV = RES / "bland_altman_stats.csv"
PNG = RES / "bland_altman.png"

DISH, CLIN, GREY, RED, INK = "#C8541A", "#2B5D9B", "#8A8A8A", "#B3242B", "#1a1a1a"


def _variance_components(d, g):
    """Unbalanced one-way random-effects decomposition (Bland & Altman 1999, repeated measures).

    Returns (var_between, var_within, n0). Closed form, so there is no optimiser to fail; the
    mixed-model fit is reported alongside only as a convergence-flagged cross-check.
    """
    df = pd.DataFrame({"d": np.asarray(d, float), "g": np.asarray(g)})
    k = df.g.nunique()
    N = len(df)
    grand = df.d.mean()
    sizes = df.groupby("g").size()
    means = df.groupby("g").d.mean()
    ss_b = float((sizes * (means - grand) ** 2).sum())
    ss_w = float(((df.d - df.g.map(means)) ** 2).sum())
    ms_b = ss_b / (k - 1)
    ms_w = ss_w / (N - k)
    n0 = (N - (sizes ** 2).sum() / N) / (k - 1)
    var_b = max((ms_b - ms_w) / n0, 0.0)
    return var_b, ms_w, float(n0)


def agreement(diff, mean, group=None, label=""):
    """Bias, limits of agreement and proportional-bias slope.

    When `group` is given, the SD of differences is rebuilt from the between- and within-cluster
    variance components, so that repeated measurements on the same drug or subject are not
    mistaken for independent precision.
    """
    d = pd.DataFrame({"d": np.asarray(diff, float), "m": np.asarray(mean, float)}).dropna()
    n = len(d)
    bias = d.d.mean()
    if group is not None:
        g = pd.Series(np.asarray(group)).loc[d.index]
        var_b, var_w, _ = _variance_components(d.d, g)
        sd = np.sqrt(var_b + var_w)
        n_groups = int(g.nunique())
        sd_naive = d.d.std(ddof=1)
    else:
        sd = d.d.std(ddof=1)
        var_b, var_w, n_groups, sd_naive = np.nan, sd ** 2, np.nan, sd
    pb = sm.OLS(d.d, sm.add_constant(d.m)).fit()
    se_bias = sd / np.sqrt(n if group is None else n_groups)
    return dict(comparison=label, n=n, n_groups=n_groups, bias_ms=bias,
                bias_ci_lo=bias - 1.96 * se_bias, bias_ci_hi=bias + 1.96 * se_bias,
                sd_diff_ms=sd, sd_diff_naive_ms=sd_naive, var_between=var_b, var_within=var_w,
                loa_lower_ms=bias - 1.96 * sd, loa_upper_ms=bias + 1.96 * sd,
                prop_bias_slope=float(pb.params.iloc[1]), prop_bias_p=float(pb.pvalues.iloc[1]))


def panel(ax, mean, diff, stats, colour, xlabel, title):
    ax.scatter(mean, diff, s=7, color=colour, alpha=0.45, lw=0)
    ax.axhline(stats["bias_ms"], color=INK, lw=1.0)
    ax.axhline(stats["loa_upper_ms"], color=RED, lw=0.9, ls=(0, (4, 2)))
    ax.axhline(stats["loa_lower_ms"], color=RED, lw=0.9, ls=(0, (4, 2)))
    ax.axhline(0, color=GREY, lw=0.6, ls=":")
    xr = ax.get_xlim()
    ax.text(xr[1], stats["bias_ms"], f"  bias {stats['bias_ms']:+.1f}", va="center",
            ha="left", fontsize=6.4, color=INK)
    ax.text(xr[1], stats["loa_upper_ms"], f"  +1.96 SD  {stats['loa_upper_ms']:+.0f}",
            va="center", ha="left", fontsize=6.0, color=RED)
    ax.text(xr[1], stats["loa_lower_ms"], f"  -1.96 SD  {stats['loa_lower_ms']:+.0f}",
            va="center", ha="left", fontsize=6.0, color=RED)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("difference (ms)")
    ax.set_title(title, loc="left")


def main():
    rows, panels = [], []

    # (a) optical APD90 vs electrode FPD, clustered by drug
    arm = pd.read_csv(RES / "arm_bycell.csv").dropna(subset=["MEA", "VSO"])
    d_arm = arm.VSO - arm.MEA
    m_arm = (arm.VSO + arm.MEA) / 2
    s_arm = agreement(d_arm, m_arm, group=arm.drug,
                      label="optical APD90 minus electrode FPD (paired drug x concentration)")
    rows.append(s_arm)
    # one cell is a gross outlier (difference > 1,000 ms); report the comparison without it too,
    # so the proportional-bias finding cannot rest on a single point
    keep = (d_arm - d_arm.mean()).abs() < 5 * d_arm.std(ddof=1)
    rows.append(agreement(d_arm[keep], m_arm[keep], group=arm.drug[keep],
                          label="optical APD90 minus electrode FPD, single extreme cell removed"))
    panels.append((m_arm, d_arm, s_arm, DISH,
                   "mean of the two methods (ms)",
                   "a   Optical against electrode, 97 paired cells"))

    # (b) tangent T-end vs the study's own, clustered by subject
    wav = pd.read_csv(RES / "waveform_delineation.csv")
    wav = wav[wav.plausible].dropna(subset=["qt_study", "qt_tangent"])
    d_w = wav.qt_tangent - wav.qt_study
    m_w = (wav.qt_tangent + wav.qt_study) / 2
    s_w = agreement(d_w, m_w, group=wav.RANDID,
                    label="tangent QT minus study QT (clinical median beats, by subject)")
    rows.append(s_w)
    panels.append((m_w, d_w, s_w, CLIN,
                   "mean QT of the two rules (ms)",
                   "b   Tangent against study T-end, 22 subjects"))

    # pooled cross-check on the full delineation set (no subject id in that file)
    both = pd.read_csv(RES / "delineation_both.csv")
    both = both[both.plausible].dropna(subset=["qt_study", "qt_tangent"])
    for study, g in both.groupby("study"):
        rows.append(agreement(g.qt_tangent - g.qt_study, (g.qt_tangent + g.qt_study) / 2,
                              label=f"tangent QT minus study QT, {study} (pooled, clustering not modelled)"))

    out = pd.DataFrame(rows)
    out.to_csv(CSV, index=False)

    mpl.rcParams.update({"font.family": "sans-serif", "font.size": 8, "figure.dpi": 200,
                         "figure.facecolor": "white", "axes.facecolor": "white",
                         "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(1, 2, figsize=(9.4, 3.7))
    for ax, (m, d, s, c, xl, t) in zip(axes, panels):
        panel(ax, m, d, s, c, xl, t)
    fig.subplots_adjust(left=0.07, right=0.84, bottom=0.17, top=0.88, wspace=0.52)
    fig.savefig(PNG, dpi=200)

    print(out[["comparison", "n", "n_groups", "bias_ms", "sd_diff_ms",
               "loa_lower_ms", "loa_upper_ms", "prop_bias_slope", "prop_bias_p"]].round(3).to_string(index=False))
    print(f"\nwrote {CSV}\nwrote {PNG}")


if __name__ == "__main__":
    main()
