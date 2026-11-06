"""Figure 2 of docs/PROJECT_REPORT.md - the dish2ecg processing pipeline, drawn.

One correction cascade is applied to in-vitro, clinical and in-silico data. Each stage is
illustrated with this repository's own numbers where raw values exist, and with a labelled
schematic where only summaries do (the recorded-signal row). Writes results/pipeline_diagram.png.

    python scripts/make_pipeline_diagram.py
"""
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
OUT = RES / "pipeline_diagram.png"

DISH, CLIN, SIL, RED, GREY, INK = "#C8541A", "#2B5D9B", "#4A4A4A", "#B3242B", "#8A8A8A", "#1a1a1a"
COL = {"dish": DISH, "clin": CLIN, "sil": SIL}
COLS = {"dish": (13.5, 41.5), "clin": (43.5, 71.5), "sil": (73.5, 99.0)}
HEAD = {"dish": "In vitro · the dish", "clin": "In vivo · the clinic", "sil": "In silico"}
L, R, B, T = 0.005, 0.995, 0.005, 0.995

ROWS = [86.5, 69.0, 52.5, 34.0]
HGT = [15.0, 13.0, 12.0, 16.0]
STEP = ["Recorded signal\nand endpoint", "Rate\ncorrection",
        "Baseline, then\nplacebo  (Δ, ΔΔ)", "Concentration–\nresponse model"]


# ---------------------------------------------------------------- schematics
def field_potential(t):
    """Extracellular field potential: fast depolarisation spike, slow repolarisation wave."""
    spike = -0.95 * np.exp(-((t - 8) ** 2) / 6) + 0.38 * np.exp(-((t - 16) ** 2) / 14)
    wave = 0.30 * np.exp(-((t - 350) ** 2) / (2 * 55 ** 2))
    return spike + wave - 0.02 * np.exp(-((t - 120) ** 2) / (2 * 90 ** 2))


def ecg_beat(t):
    p = 0.11 * np.exp(-((t - 150) ** 2) / (2 * 22 ** 2))
    q = -0.09 * np.exp(-((t - 292) ** 2) / (2 * 7 ** 2))
    r = 1.35 * np.exp(-((t - 305) ** 2) / (2 * 8 ** 2))
    s = -0.16 * np.exp(-((t - 322) ** 2) / (2 * 9 ** 2))
    tw = 0.32 * np.exp(-((t - 600) ** 2) / (2 * 55 ** 2))
    return p + q + r + s + tw


def action_potential(t):
    up = 1 / (1 + np.exp(-(t - 4) / 0.8))
    plateau = 1 - 0.30 * (1 / (1 + np.exp(-(t - 90) / 55)))
    repol = 1 / (1 + np.exp((t - 245) / 22))
    return -72 + 104 * up * np.maximum(plateau * repol, 0)


# ---------------------------------------------------------------- helpers
def rect(fig, x0, x1, yc, h):
    fx = lambda x: L + (x / 100) * (R - L)
    fy = lambda y: B + (y / 100) * (T - B)
    return fig.add_axes([fx(x0), fy(yc - h / 2), fx(x1) - fx(x0), fy(yc + h / 2) - fy(yc - h / 2)])


def panel(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_linewidth(0.6)
        ax.spines[s].set_color("#999999")
    ax.tick_params(labelsize=5.2, length=1.8, width=0.5, pad=1.2, colors="#555555")
    return ax


def trace_panel(ax, c):
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    return ax


def box(ax, x0, x1, yc, h, ec, fc="white", ls="-", lw=1.0):
    ax.add_patch(FancyBboxPatch((x0, yc - h / 2), x1 - x0, h, boxstyle="round,pad=0.4,rounding_size=1.2",
                                linewidth=lw, edgecolor=ec, facecolor=fc, linestyle=ls, zorder=1))


def arrow(ax, x, y0, y1, c):
    ax.add_patch(FancyArrowPatch((x, y0), (x, y1), arrowstyle="-|>", mutation_scale=9,
                                 linewidth=0.9, color=c, zorder=2))


def bracket(ax, x0, x1, y, label, c, fs=5.6, dy=0.0):
    ax.annotate("", xy=(x0, y), xytext=(x1, y),
                arrowprops=dict(arrowstyle="<|-|>", color=c, lw=0.9, mutation_scale=6, shrinkA=0, shrinkB=0))
    ax.text((x0 + x1) / 2, y + dy, label, ha="center", va="bottom", fontsize=fs, color=c)


def cap(bg, x0, x1, y, txt, c=INK, fs=5.9):
    bg.text((x0 + x1) / 2, y, txt, ha="center", va="center", fontsize=fs, color=c, linespacing=1.4)


# ---------------------------------------------------------------- data
def load():
    d = {}
    d["millard"] = pd.read_csv(RES / "invitro_correction_validation.csv")
    cc = pd.read_csv(RES / "clinical_ddqtc.csv")
    e = cc[cc.study == "ecgrdvq"]
    d["free"] = e[(e.arm == "Placebo") | (e.is_baseline)][["RR", "QT"]].dropna()
    dof = e[e.arm == "Dofetilide"]
    d["dof"] = dof[(~dof.is_baseline)][["conc_ngml", "ddQTcF"]].dropna()
    tpt = dof[~dof.is_baseline].groupby("TPT").ddQTcF.mean().idxmax()
    on = dof[(~dof.is_baseline) & (dof.TPT == tpt)]
    base = dof[dof.is_baseline].groupby("RANDID").QTcF.mean()
    plac = e[e.arm == "Placebo"]
    pb = plac[plac.is_baseline].groupby("RANDID").QTcF.mean()
    po = plac[(~plac.is_baseline) & (plac.TPT == tpt)]
    d["wf"] = dict(base=float(base.reindex(on.RANDID).mean()), on=float(on.QTcF.mean()),
                   delta=float((on.set_index("RANDID").QTcF - base).dropna().mean()),
                   plac=float((po.set_index("RANDID").QTcF - pb).dropna().mean()),
                   dd=float(on.ddQTcF.mean()), n=int(len(on)), tpt=float(tpt))
    iv = pd.read_csv(RES / "invitro_concentration_response.csv")
    r = iv[(iv.drug.str.lower() == "dofetilide") & (iv.analysis == "all_concentrations") &
           (iv.endpoint == "ddFPDc")].iloc[0]
    d["dishfit"] = dict(slope=float(r.slope), intercept=float(r.intercept),
                        slope_se=float(r.slope_se), intercept_se=float(r.intercept_se),
                        centre=float(r.at_concentration), n=int(r.n_obs), g=int(r.n_groups))
    cq = pd.read_csv(RES / "concentration_qtc_results.csv")
    q = cq[(cq.endpoint == "ddQTcF") & (cq.study == "ecgrdvq") & (cq.arm == "Dofetilide")].iloc[0]
    d["clinfit"] = dict(slope=float(q.slope), intercept=float(q.intercept),
                        slope_se=float(q.slope_se), intercept_se=float(q.intercept_se),
                        centre=float(q.at_concentration), g=int(q.n_groups),
                        ci_hi=float(q.effect_ci90_upper), exceeds=bool(q.exceeds_10ms))
    od = pd.read_csv(RES / "ord_apd90.csv")
    o = od[(od.drug == "dofetilide") & (od.cell == "mid")].sort_values("multiple")
    # Beyond ~3.8x the model stops repolarising: apd90 is NaN or the "resting" potential is
    # depolarised and the detector returns a meaningless value. The `ok` flag does not catch
    # this, so filter on the physiology instead.
    d["ord"] = o[o.apd90.notna() & (o.v_rest < -60)]
    allcells = od[(od.drug == "dofetilide") & od.apd90.notna() & (od.v_rest < -60)]
    d["ord_all"] = allcells
    d["ord_fail"] = float(o[o.apd90.isna() | (o.v_rest >= -60)].multiple.min())
    d["bridge"] = pd.read_csv(RES / "dish_vs_clinic_bridge.csv")
    return d


# ---------------------------------------------------------------- stage rows
def row_signal(fig, bg, d):
    y, h = ROWS[0], HGT[0]
    for k in COLS:
        box(bg, *COLS[k], y, h, COL[k])
    ph = h - 6.2
    py = y + 1.1

    x0, x1 = COLS["dish"]
    ax = trace_panel(rect(fig, x0 + 2.2, x1 - 2.2, py, ph), DISH)
    t = np.linspace(-40, 560, 1400)
    ax.plot(t, field_potential(t), color=DISH, lw=1.0)
    ax.axhline(0, color="#CCCCCC", lw=0.5, zorder=0)
    ax.set_ylim(-1.25, 0.72); ax.set_xlim(-40, 560)
    bracket(ax, 8, 350, 0.50, "FPD", DISH)
    ax.text(8, -1.12, "depolarisation\nspike", fontsize=4.9, color=GREY, ha="left", va="bottom", linespacing=1.3)
    ax.text(455, -0.14, "repolarisation\nwave peak", fontsize=4.9, color=GREY, ha="left", va="top", linespacing=1.3)
    cap(bg, x0, x1, y - h / 2 + 2.3, "field potential under a monolayer  ·  schematic\n"
        "28 drugs · 10 sites · 15 datasets · 8,859 wells")

    x0, x1 = COLS["clin"]
    ax = trace_panel(rect(fig, x0 + 2.2, x1 - 2.2, py, ph), CLIN)
    t = np.linspace(0, 900, 1800)
    ax.plot(t, ecg_beat(t), color=CLIN, lw=1.0)
    ax.axhline(0, color="#CCCCCC", lw=0.5, zorder=0)
    ax.set_ylim(-0.55, 1.95); ax.set_xlim(0, 900)
    bracket(ax, 288, 690, 1.45, "QT", CLIN)
    ax.text(150, 0.22, "P", fontsize=5.2, color=GREY, ha="center")
    ax.text(305, 1.44, "R", fontsize=5.2, color=GREY, ha="center")
    ax.text(600, 0.42, "T", fontsize=5.2, color=GREY, ha="center")
    cap(bg, x0, x1, y - h / 2 + 2.3, "12-lead ECG, median beat  ·  schematic\n"
        "104 volunteers · 11 molecules · 3 PhysioNet studies")

    x0, x1 = COLS["sil"]
    ax = trace_panel(rect(fig, x0 + 2.2, x1 - 2.2, py, ph), SIL)
    t = np.linspace(-20, 420, 1400)
    v = action_potential(t)
    ax.plot(t, v, color=SIL, lw=1.0)
    thr = -72 + 0.1 * (v.max() - (-72))
    ax.axhline(thr, color="#CCCCCC", lw=0.5, ls=(0, (3, 2)), zorder=0)
    ax.set_ylim(-88, 60); ax.set_xlim(-20, 420)
    i90 = np.where(v <= thr)[0]
    t90 = t[i90[i90 > 400][0]] if len(i90[i90 > 400]) else 265
    bracket(ax, 4, t90, 42, "APD$_{90}$", SIL)
    ax.text(300, thr + 3, "90 % repolarised", fontsize=4.9, color=GREY, ha="left", va="bottom")
    cap(bg, x0, x1, y - h / 2 + 2.3, "simulated ventricular action potential\n"
        "28 drugs × 7 currents · O'Hara–Rudy / ORd-CiPA")

    # guardrail, drawn where it fires
    gx0, gx1 = COLS["dish"]
    bg.text((gx0 + gx1) / 2, y + h / 2 - 1.6, "✗  expressing this as QT in ms raises DomainError",
            ha="center", va="center", fontsize=5.8, color=RED)


def row_rate(fig, bg, d):
    y, h = ROWS[1], HGT[1]
    for k in COLS:
        box(bg, *COLS[k], y, h, COL[k])
    ph, py = h - 6.0, y + 0.9

    m = d["millard"]
    x0, x1 = COLS["dish"]
    ax = panel(rect(fig, x0 + 4.2, x1 - 2.2, py, ph))
    for _, r in m.iterrows():
        ax.plot([r.BP_ms, r.BP_ms], [r.FPD_ms, r.FPDc_dish2ecg], color="#D9B49C", lw=0.6, zorder=1)
    ax.scatter(m.BP_ms, m.FPD_ms, s=9, facecolor="none", edgecolor=DISH, lw=0.7, zorder=3, label="FPD")
    ax.scatter(m.BP_ms, m.FPDc_dish2ecg, s=9, color=DISH, lw=0, zorder=3, label="FPDc")
    ax.legend(fontsize=4.8, frameon=False, loc="upper left", handlelength=0.7,
              borderpad=0.1, handletextpad=0.3, labelspacing=0.2)
    ax.set_xlabel("beat period (ms)", fontsize=5.4, labelpad=1.2)
    ax.set_ylabel("FPD (ms)", fontsize=5.4, labelpad=1.2)
    cap(bg, x0, x1, y - h / 2 + 2.1, "FPDc = FPD / BP$^{1/3}$   ·   each of 18 studies, Millard 2018")

    f = d["free"]
    x0, x1 = COLS["clin"]
    ax = panel(rect(fig, x0 + 4.2, x1 - 2.2, py, ph))
    ax.scatter(f.RR, f.QT, s=2.2, color=CLIN, alpha=0.30, lw=0)
    rr = np.linspace(f.RR.min(), f.RR.max(), 100)
    qtcf_med = float((f.QT / (f.RR / 1000) ** (1 / 3)).median())
    ax.plot(rr, qtcf_med * (rr / 1000) ** (1 / 3), color=GREY, lw=0.9, ls=(0, (3, 2)))
    ax.text(0.96, 0.06, "dashed: RR$^{1/3}$ through\nthe median QTcF", transform=ax.transAxes,
            ha="right", va="bottom", fontsize=4.6, color=GREY, linespacing=1.3)
    ax.set_xlabel("RR (ms)", fontsize=5.4, labelpad=1.2)
    ax.set_ylabel("QT (ms)", fontsize=5.4, labelpad=1.2)
    cap(bg, x0, x1, y - h / 2 + 2.1,
        f"QTcF = QT / RR$^{{1/3}}$   ·   {len(f):,} drug-free records")

    x0, x1 = COLS["sil"]
    box(bg, x0 + 2.0, x1 - 2.0, y + 0.6, h - 5.4, SIL, "#F5F5F5", ls=(0, (3, 2)), lw=0.8)
    bg.text((x0 + x1) / 2, y + 0.6, "not applicable\n\neach cell is paced at a\nfixed cycle length, so there\nis no rate to correct for",
            ha="center", va="center", fontsize=5.9, color="#555555", linespacing=1.5)
    cap(bg, x0, x1, y - h / 2 + 2.1, "1,000 ms pacing, pre-paced to steady state")


def row_subtract(fig, bg, d):
    y, h = ROWS[2], HGT[2]
    for k in COLS:
        box(bg, *COLS[k], y, h, COL[k])
    ph, py = h - 5.6, y + 0.9
    w = d["wf"]

    x0, x1 = COLS["clin"]
    ax = panel(rect(fig, x0 + 4.2, x1 - 2.2, py, ph))
    vals = [w["base"], w["on"], w["delta"], w["plac"], w["dd"]]
    lab = ["baseline", "on drug", "Δ", "placebo Δ", "ΔΔ"]
    cols = ["#BBBBBB", CLIN, CLIN, "#BBBBBB", CLIN]
    ax.bar(range(5), vals, color=cols, width=0.62, lw=0)
    for i, v in enumerate(vals):
        ax.text(i, v + (5 if v >= 0 else -5), f"{v:+.0f}" if i >= 2 else f"{v:.0f}",
                ha="center", va="bottom" if v >= 0 else "top", fontsize=5.0, color=INK)
    ax.axhline(0, color="#999999", lw=0.5)
    ax.set_xticks(range(5)); ax.set_xticklabels(lab, fontsize=4.9)
    ax.set_ylabel("QTcF (ms)", fontsize=5.4, labelpad=1.2)
    ax.set_ylim(-40, 560)
    cap(bg, x0, x1, y - h / 2 + 2.0,
        f"dofetilide, {w['n']} subjects at {w['tpt']:g} h   ·   ΔΔ = Δ − placebo Δ")

    x0, x1 = COLS["dish"]
    box(bg, x0 + 2.0, x1 - 2.0, y + 0.6, h - 5.0, DISH, "#FBF2EC", ls=(0, (3, 2)), lw=0.8)
    bg.text((x0 + x1) / 2, y + 0.6,
            "already applied in the source\n\nthe published workbook ships ΔΔFPDc:\n"
            "baseline-, vehicle- and rate-corrected\nby the study, verified here, not re-applied",
            ha="center", va="center", fontsize=5.9, color="#6B4428", linespacing=1.5)
    cap(bg, x0, x1, y - h / 2 + 2.0, "verified against 420 published panels: r = 0.999993")

    x0, x1 = COLS["sil"]
    box(bg, x0 + 2.0, x1 - 2.0, y + 0.6, h - 5.0, SIL, "#F5F5F5", ls=(0, (3, 2)), lw=0.8)
    bg.text((x0 + x1) / 2, y + 0.6,
            "not applicable\n\nthe drug-free run of the same\nmodel is the control, so there is\nno baseline or placebo to subtract",
            ha="center", va="center", fontsize=5.9, color="#555555", linespacing=1.5)
    cap(bg, x0, x1, y - h / 2 + 2.0, "ΔAPD$_{90}$ is taken against the same cell, drug-free")


def row_model(fig, bg, d):
    y, h = ROWS[3], HGT[3]
    for k in COLS:
        box(bg, *COLS[k], y, h, COL[k])
    ph, py = h - 7.4, y + 1.9

    fit = d["dishfit"]
    x0, x1 = COLS["dish"]
    ax = panel(rect(fig, x0 + 4.4, x1 - 2.2, py, ph))
    lg = np.linspace(-1.0, 0.75, 100)
    yhat = fit["intercept"] + fit["slope"] * (lg - fit["centre"])
    tcrit = stats.t.ppf(0.95, fit["g"] - 1)          # two-sided 90%, as the cascade reports
    se = np.sqrt(fit["intercept_se"] ** 2 + ((lg - fit["centre"]) * fit["slope_se"]) ** 2)
    ax.fill_between(lg, yhat - tcrit * se, yhat + tcrit * se, color=DISH, alpha=0.16, lw=0)
    ax.plot(lg, yhat, color=DISH, lw=1.2)
    ax.axhline(0, color="#CCCCCC", lw=0.5)
    ax.axvline(np.log10(1.074), color=GREY, lw=0.6, ls=(0, (2, 2)))
    ax.scatter([np.log10(1.074)], [132.0], s=16, color=DISH, zorder=5, edgecolor="white", lw=0.5)
    ax.annotate("+132 ms", xy=(np.log10(1.074), 132.0), xytext=(-0.92, 182),
                fontsize=5.2, color=DISH,
                arrowprops=dict(arrowstyle="-", color=DISH, lw=0.6, shrinkA=0, shrinkB=2))
    ax.set_xlabel("log$_{10}$(× free C$_{max}$)", fontsize=5.4, labelpad=1.2)
    ax.set_ylabel("ΔΔFPDc (ms)", fontsize=5.4, labelpad=1.2)
    ax.set_ylim(-60, 260)
    ax.text(0.97, 0.05, "no in-vitro equivalent of\nthe 10 ms rule exists", transform=ax.transAxes,
            ha="right", va="bottom", fontsize=4.7, color=GREY, linespacing=1.3)
    cap(bg, x0, x1, y - h / 2 + 1.5,
        f"mixed model · 90 % band · {fit['n']} wells, {fit['g']} datasets")

    dof = d["dof"]
    x0, x1 = COLS["clin"]
    ax = panel(rect(fig, x0 + 4.4, x1 - 2.2, py, ph))
    ax.scatter(dof.conc_ngml, dof.ddQTcF, s=2.4, color=CLIN, alpha=0.30, lw=0)
    cf = d["clinfit"]
    xs = np.linspace(0, dof.conc_ngml.max(), 100)
    yhat = cf["intercept"] + cf["slope"] * (xs - cf["centre"])
    tcrit = stats.t.ppf(0.95, cf["g"] - 1)
    se = np.sqrt(cf["intercept_se"] ** 2 + ((xs - cf["centre"]) * cf["slope_se"]) ** 2)
    ax.fill_between(xs, yhat - tcrit * se, yhat + tcrit * se, color=CLIN, alpha=0.18, lw=0)
    ax.plot(xs, yhat, color=CLIN, lw=1.2)
    ax.axhline(0, color="#CCCCCC", lw=0.5)
    ax.axhline(10, color=RED, lw=0.8, ls=(0, (4, 2)))
    ax.text(0.985, 10, "10 ms, ICH E14 ", transform=ax.get_yaxis_transform(),
            ha="right", va="bottom", fontsize=4.9, color=RED)
    ax.scatter([2.70], [73.1], s=16, color=CLIN, zorder=5, edgecolor="white", lw=0.5)
    ax.annotate(f"+73 ms\nupper 90 % bound {cf['ci_hi']:.1f} ms\n{'above' if cf['exceeds'] else 'below'} the 10 ms rule",
                xy=(2.70, 73.1), xytext=(0.10, 128), fontsize=5.0, color=CLIN, linespacing=1.35,
                arrowprops=dict(arrowstyle="-", color=CLIN, lw=0.6, shrinkA=0, shrinkB=2))
    ax.set_xlabel("plasma concentration (ng/mL)", fontsize=5.4, labelpad=1.2)
    ax.set_ylabel("ΔΔQTcF (ms)", fontsize=5.4, labelpad=1.2)
    cap(bg, x0, x1, y - h / 2 + 1.5,
        f"mixed model · 90 % band · {len(dof)} records, 22 subjects")

    o = d["ord"]
    x0, x1 = COLS["sil"]
    ax = panel(rect(fig, x0 + 4.4, x1 - 2.2, py, ph))
    spread = d["ord_all"].pivot_table(index="multiple", columns="cell", values="dAPD90").dropna()
    ax.fill_between(spread.index, spread.min(axis=1), spread.max(axis=1),
                    color=SIL, alpha=0.16, lw=0)
    ax.plot(o.multiple, o.dAPD90, color=SIL, lw=1.2, marker="o", ms=2.0, mfc=SIL, mec="none")
    ax.axhline(0, color="#CCCCCC", lw=0.5)
    ax.axvline(d["ord_fail"], color=RED, lw=0.7, ls=(0, (2, 2)))
    ax.text(d["ord_fail"] * 0.86, ax.get_ylim()[1] * 0.14,
            "beyond here the model\nstops repolarising", fontsize=4.8, color=RED,
            ha="right", va="bottom", linespacing=1.3)
    ax.set_xscale("log")
    ax.set_xlim(0.025, d["ord_fail"] * 3.2)
    ax.set_xlabel("× free C$_{max}$", fontsize=5.4, labelpad=1.2)
    ax.set_ylabel("ΔAPD$_{90}$ (ms)", fontsize=5.4, labelpad=1.2)
    cap(bg, x0, x1, y - h / 2 + 1.5, "mid-myocardial cell · band = range over the 3 cell types")


def row_tail(fig, bg, d):
    # shared exposure scale
    box(bg, 13.5, 99.0, 21.8, 5.4, INK, "#F0F0F0")
    bg.text(56.25, 22.9, "Shared exposure scale:  multiples of the geometric-mean $observed$ clinical free C$_{max}$",
            ha="center", va="center", fontsize=6.6, color=INK)
    bg.text(56.25, 20.6, "ng/mL  →  nM (molecular weight)  →  × fraction unbound  →  ÷ free C$_{max}$",
            ha="center", va="center", fontsize=6.2, color="#555555")
    for k in COLS:
        arrow(bg, sum(COLS[k]) / 2, ROWS[3] - HGT[3] / 2, 24.6, COL[k])

    # comparison
    y, h = 9.6, 12.8
    box(bg, 13.5, 99.0, y, h, INK)
    arrow(bg, 56.25, 19.1, y + h / 2 + 0.4, INK)
    br = d["bridge"].copy()
    cq = pd.read_csv(RES / "concentration_qtc_results.csv")
    cq = cq[(cq.endpoint == "ddQTcF") & (cq.study == "ecgrdvq")]
    e14 = {r.arm.split()[0].lower(): bool(r.exceeds_10ms) for _, r in cq.iterrows()}
    br["e14"] = [e14.get(n.lower(), None) for n in br.drug]
    ax = panel(rect(fig, 20.0, 62.0, y - 0.4, h - 5.6))
    idx = np.arange(len(br))
    ax.axhline(0, color="#BBBBBB", lw=0.6)
    ax.axhline(10, color=RED, lw=0.8, ls=(0, (4, 2)), zorder=2)
    ax.text(len(br) - 0.55, 12, "10 ms, ICH E14", fontsize=5.0, color=RED, ha="right", va="bottom")
    ax.scatter(idx - 0.13, br.clin_ddQTcF_ms, s=20, color=CLIN, zorder=4, label="clinic  ΔΔQTcF")
    ax.scatter(idx + 0.13, br.dish_ddFPDc_ms, s=20, color=DISH, marker="s", zorder=4, label="dish  ΔΔFPDc")
    for i, r in br.iterrows():
        ax.plot([i - 0.13, i + 0.13], [r.clin_ddQTcF_ms, r.dish_ddFPDc_ms],
                color="#CCCCCC", lw=0.8, zorder=1)
    ax.set_xticks(idx)
    ax.set_xticklabels([n[:4].lower() + "." for n in br.drug], fontsize=5.2)
    ax.set_xlim(-0.5, len(br) - 0.5)
    lo = min(br.clin_ddQTcF_ms.min(), br.dish_ddFPDc_ms.min())
    hi = max(br.clin_ddQTcF_ms.max(), br.dish_ddFPDc_ms.max())
    ax.set_ylim(lo - 0.28 * (hi - lo), hi + 0.22 * (hi - lo))
    ax.set_ylabel("change (ms)", fontsize=5.4, labelpad=1.2)
    ax.legend(fontsize=5.0, frameon=False, loc="upper right", handlelength=0.9,
              borderpad=0.1, handletextpad=0.4, labelspacing=0.25)
    ax.annotate("verapamil inverts", xy=(3.13, -86.3), xytext=(2.05, -62), fontsize=5.2, color=RED,
                ha="right", va="center",
                arrowprops=dict(arrowstyle="->", color=RED, lw=0.7, shrinkA=0, shrinkB=3))

    bg.text(65.0, y + 5.0, "Comparison at matched exposure", ha="left", va="center",
            fontsize=7.4, color=INK, weight="bold")
    bg.text(65.0, y + 3.4,
            "Four drugs are the entire overlap between the public\n"
            "in-vitro panel and the public clinical corpora.\n\n"
            "Sign agreement:  dish 3 of 4,  in silico 4 of 4.\n"
            "On the three it gets right, the dish over-reads 1.8–3.5×.\n\n"
            "ICH E14 decision, upper 90 % bound above 10 ms:\n"
            "dofetilide 77.7, quinidine 80.5, ranolazine 14.8 — yes.\n"
            "Verapamil 7.0 — no. The dish has no matching rule.",
            ha="left", va="top", fontsize=6.0, color="#333333", linespacing=1.6)


def main():
    mpl.rcParams.update({"font.family": "sans-serif", "font.size": 8.0, "figure.dpi": 200,
                         "figure.facecolor": "white", "axes.facecolor": "white",
                         "mathtext.fontset": "dejavusans", "mathtext.default": "regular"})
    d = load()
    fig = plt.figure(figsize=(11.0, 14.2))
    bg = fig.add_axes([L, B, R - L, T - B])
    bg.set_xlim(0, 100); bg.set_ylim(0, 100); bg.axis("off")
    bg.patch.set_alpha(0)

    for k, (x0, x1) in COLS.items():
        bg.text((x0 + x1) / 2, 97.3, HEAD[k], ha="center", va="center",
                fontsize=9.2, weight="bold", color=COL[k])
        bg.plot([x0, x1], [95.4, 95.4], color=COL[k], lw=1.6, solid_capstyle="butt")
    bg.text(6.5, 97.3, "Shared\ncascade step", ha="center", va="center",
            fontsize=7.2, weight="bold", color=INK, linespacing=1.4)
    for yy, s in zip(ROWS, STEP):
        bg.text(6.5, yy, s, ha="center", va="center", fontsize=6.8, color="#333333", linespacing=1.5)

    for i in range(len(ROWS) - 1):
        for k in ("dish", "clin"):
            arrow(bg, sum(COLS[k]) / 2, ROWS[i] - HGT[i] / 2, ROWS[i + 1] + HGT[i + 1] / 2, COL[k])
    arrow(bg, sum(COLS["sil"]) / 2, ROWS[0] - HGT[0] / 2, ROWS[1] + HGT[1] / 2, SIL)
    arrow(bg, sum(COLS["sil"]) / 2, ROWS[1] - HGT[1] / 2, ROWS[2] + HGT[2] / 2, SIL)
    arrow(bg, sum(COLS["sil"]) / 2, ROWS[2] - HGT[2] / 2, ROWS[3] + HGT[3] / 2, SIL)

    row_signal(fig, bg, d); row_rate(fig, bg, d); row_subtract(fig, bg, d)
    row_model(fig, bg, d); row_tail(fig, bg, d)

    bg.text(13.5, 1.6,
            "The dish axis is nominal concentration; the clinical axis is measured free concentration. They are not the same quantity.",
            ha="left", va="center", fontsize=6.0, color=RED, style="italic")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=200)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
