"""Were the processing choices right? Each one, against its alternative.

Every number this project reports depends on four choices that could have gone
another way: the rate-correction exponent, the shape of the concentration
predictor, the random-effect structure, and whether the placebo adjustment is
pre-formed or estimated inside the model. This script refits the headline
quantity under each alternative and reports how far it moves.

The point is not to find the "true" number. It is to show which conclusions are
choice-dependent and which are not.

Writes results/methodology_sensitivity.csv.

    python scripts/methodology_sensitivity.py
"""
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
OUT = RES / "methodology_sensitivity.csv"

ARMS = {"Dofetilide": "Dofetilide", "Quinidine": "Quinidine Sulph",
        "Ranolazine": "Ranolazine", "Verapamil": "Verapamil HCL"}


def _fit_at(df, ycol, xcol, group, at, re_formula=None):
    """Mixed model of y on x, centred at `at`; returns effect there and its 90 % CI."""
    d = df[[ycol, xcol, group]].dropna().rename(columns={ycol: "_y", xcol: "_x", group: "_g"})
    d["_xc"] = d._x - at
    scale = max(float(np.std(d._xc)), 1e-12)
    d["_xs"] = d._xc / scale
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fit = smf.mixedlm("_y ~ _xs", data=d, groups="_g", re_formula=re_formula).fit(reml=True)
    from scipy import stats
    t = stats.t.ppf(0.95, d._g.nunique() - 1)
    b, se = float(fit.params["Intercept"]), float(fit.bse["Intercept"])
    return b, b - t * se, b + t * se, int(len(d)), int(d._g.nunique())


def clinical(rows):
    clin = pd.read_csv(RES / "clinical_ddqtc.csv")
    e = clin[clin.study == "ecgrdvq"]
    cq = pd.read_csv(RES / "concentration_qtc_results.csv")
    for drug, arm in ARMS.items():
        ref = cq[(cq.endpoint == "ddQTcF") & (cq.study == "ecgrdvq") & (cq.arm == arm)].iloc[0]
        cmax = float(ref.at_concentration)
        on = e[(e.arm == arm) & (~e.is_baseline)]

        # 1. rate-correction choice
        for endpoint, label in (("ddQTcF", "Fridericia (used)"),
                                ("ddQTcB", "Bazett"),
                                ("ddQTcI", "individual exponent")):
            if endpoint not in on or on[endpoint].notna().sum() < 20:
                continue
            b, lo, hi, n, g = _fit_at(on, endpoint, "conc_ngml", "RANDID", cmax)
            rows.append(dict(arm="in vivo", drug=drug, question="rate correction",
                             variant=label, effect_ms=b, ci_lo=lo, ci_hi=hi,
                             n=n, groups=g, exceeds_10ms=bool(hi > 10)))

        # 2. random structure
        b, lo, hi, n, g = _fit_at(on, "ddQTcF", "conc_ngml", "RANDID", cmax, re_formula="1 + _xs")
        rows.append(dict(arm="in vivo", drug=drug, question="random structure",
                         variant="intercept + slope", effect_ms=b, ci_lo=lo, ci_hi=hi,
                         n=n, groups=g, exceeds_10ms=bool(hi > 10)))

        # 3. pre-formed double delta vs white-paper-style model
        base = e[(e.arm == arm) & (e.is_baseline)].groupby("RANDID").QTcF.mean()
        pb = e[(e.arm == "Placebo") & (e.is_baseline)].groupby("RANDID").QTcF.mean()
        plac = e[(e.arm == "Placebo") & (~e.is_baseline)]
        d_on = on.assign(dQTcF=on.QTcF - on.RANDID.map(base), trt=1,
                         conc=on.conc_ngml, b0=on.RANDID.map(base))
        d_pl = plac.assign(dQTcF=plac.QTcF - plac.RANDID.map(pb), trt=0,
                           conc=0.0, b0=plac.RANDID.map(pb))
        both = pd.concat([d_on, d_pl], ignore_index=True).dropna(subset=["dQTcF", "conc", "TPT", "b0"])
        both["b0c"] = both.b0 - float(both.b0.mean())
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                wp = smf.mixedlm("dQTcF ~ conc + trt + b0c + C(TPT)", data=both,
                                 groups=both.RANDID).fit(reml=True)
            eff = float(wp.params["conc"]) * cmax + float(wp.params["trt"])
            rows.append(dict(arm="in vivo", drug=drug, question="placebo adjustment",
                             variant="white-paper model (estimated)", effect_ms=eff,
                             ci_lo=np.nan, ci_hi=np.nan, n=int(len(both)),
                             groups=int(both.RANDID.nunique()), exceeds_10ms=None))
        except Exception as exc:                                      # pragma: no cover
            rows.append(dict(arm="in vivo", drug=drug, question="placebo adjustment",
                             variant=f"white-paper model FAILED: {type(exc).__name__}",
                             effect_ms=np.nan, ci_lo=np.nan, ci_hi=np.nan,
                             n=np.nan, groups=np.nan, exceeds_10ms=None))
        rows.append(dict(arm="in vivo", drug=drug, question="placebo adjustment",
                         variant="pre-formed double delta (used)",
                         effect_ms=float(ref.effect_at_concentration),
                         ci_lo=float(ref.effect_ci90_lower), ci_hi=float(ref.effect_ci90_upper),
                         n=int(ref.n_obs), groups=int(ref.n_groups),
                         exceeds_10ms=bool(ref.exceeds_10ms)))


def invitro(rows):
    wells_path = ROOT / "data" / "cipa_myocyte_wells.csv"
    iv = pd.read_csv(RES / "invitro_concentration_response.csv")
    bridge = pd.read_csv(RES / "dish_vs_clinic_bridge.csv")
    for _, br in bridge.iterrows():
        drug, match = br.drug, float(br.match_point_folds_free_cmax)
        for analysis, label in (("all_concentrations", "all 4 concentrations (used)"),
                                ("conc_1_to_3_only", "lowest 3 concentrations")):
            sub = iv[(iv.drug == drug) & (iv.analysis == analysis) & (iv.endpoint == "ddFPDc")]
            if sub.empty:
                continue
            r = sub.iloc[0]
            from scipy import stats
            t = stats.t.ppf(0.95, int(r.n_groups) - 1)
            x = np.log10(match) - float(r.at_concentration)
            eff = float(r.intercept) + float(r.slope) * x
            se = np.sqrt(float(r.intercept_se) ** 2 + (x * float(r.slope_se)) ** 2)
            rows.append(dict(arm="in vitro", drug=drug, question="informative missingness",
                             variant=label, effect_ms=eff, ci_lo=eff - t * se, ci_hi=eff + t * se,
                             n=int(r.n_obs), groups=int(r.n_groups), exceeds_10ms=None))
        if not wells_path.exists():
            continue
        w = pd.read_csv(wells_path)
        w = w[(w.drug.str.lower() == drug.lower()) & w.ddFPDc.notna() &
              w.folds_cmax.notna() & (w.folds_cmax > 0)]
        if len(w) < 30:
            continue
        w = w.assign(logm=np.log10(w.folds_cmax),
                     dataset=w.site.astype(str) + "|" + w.cell_line.astype(str) + "|" + w.Platform.astype(str))
        # predictor shape
        b, lo, hi, n, g = _fit_at(w, "ddFPDc", "logm", "dataset", np.log10(match))
        rows.append(dict(arm="in vitro", drug=drug, question="predictor shape",
                         variant="log10 multiples (used)", effect_ms=b, ci_lo=lo, ci_hi=hi,
                         n=n, groups=g, exceeds_10ms=None))
        b, lo, hi, n, g = _fit_at(w, "ddFPDc", "folds_cmax", "dataset", match)
        rows.append(dict(arm="in vitro", drug=drug, question="predictor shape",
                         variant="linear in multiples", effect_ms=b, ci_lo=lo, ci_hi=hi,
                         n=n, groups=g, exceeds_10ms=None))
        # dataset as fixed effect instead of random
        d = w.assign(_xc=w.logm - np.log10(match))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            ols = smf.ols("ddFPDc ~ _xc + C(dataset)", data=d).fit()
        pred = d.assign(_xc=0.0)
        eff = float(ols.predict(pred).mean())
        rows.append(dict(arm="in vitro", drug=drug, question="dataset effect",
                         variant="fixed (mean over datasets)", effect_ms=eff, ci_lo=np.nan,
                         ci_hi=np.nan, n=int(len(d)), groups=int(d.dataset.nunique()),
                         exceeds_10ms=None))


def main():
    rows = []
    clinical(rows)
    invitro(rows)
    out = pd.DataFrame(rows)
    out.to_csv(OUT, index=False)
    with pd.option_context("display.width", 220, "display.max_rows", 200):
        print(out.round(2).to_string(index=False))
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
