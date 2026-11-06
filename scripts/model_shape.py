"""Is the concentration-response linear? Linear against log-linear and Emax.

The white-paper model is linear and regulators treat it as the default, but the same
guidance requires linearity to be CHECKED and an Emax model fitted when it fails. This
script does that check on the clinical arm, and states why the question does not arise
for the in-silico arm.

Structural forms, each with the same random intercept per subject, fitted by ML so the
fixed-effect structures are comparable by BIC:

    linear      f(C) = beta * C
    log-linear  f(C) = p1 * log(1 + C / p2)
    Emax        f(C) = Emax * C / (C + EC50)

p2 and EC50 are profiled over a grid; conditional on them each model is linear in its
parameters, so the same mixed-model machinery applies throughout and the comparison is
not confounded by a change of estimator.

Writes results/model_shape.csv.

    python scripts/model_shape.py
"""
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
OUT = RES / "model_shape.csv"

ARMS = {"Dofetilide": "Dofetilide", "Quinidine": "Quinidine Sulph",
        "Ranolazine": "Ranolazine", "Verapamil": "Verapamil HCL"}


def _mixed(y, x, g):
    d = pd.DataFrame({"_y": y, "_x": x, "_g": g}).dropna()
    s = max(float(np.std(d._x)), 1e-12)
    d["_xs"] = d._x / s
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fit = smf.mixedlm("_y ~ _xs", data=d, groups="_g").fit(reml=False)
    return fit, s, d


def _effect_at(fit, s, xval, df):
    b0, b1 = float(fit.params["Intercept"]), float(fit.params["_xs"])
    v = fit.cov_params()
    var = (v.loc["Intercept", "Intercept"] + (xval / s) ** 2 * v.loc["_xs", "_xs"]
           + 2 * (xval / s) * v.loc["Intercept", "_xs"])
    se = float(np.sqrt(max(var, 0.0)))
    eff = b0 + b1 * xval / s
    t = stats.t.ppf(0.95, df)
    return eff, eff - t * se, eff + t * se


def main():
    clin = pd.read_csv(RES / "clinical_ddqtc.csv")
    e = clin[clin.study == "ecgrdvq"]
    cq = pd.read_csv(RES / "concentration_qtc_results.csv")
    rows = []
    for drug, arm in ARMS.items():
        ref = cq[(cq.endpoint == "ddQTcF") & (cq.study == "ecgrdvq") & (cq.arm == arm)].iloc[0]
        cmax = float(ref.at_concentration)
        on = e[(e.arm == arm) & (~e.is_baseline)][["RANDID", "conc_ngml", "ddQTcF"]].dropna()
        y, C, g = on.ddQTcF.values, on.conc_ngml.values, on.RANDID.values
        ndf = on.RANDID.nunique() - 1
        cand = {}

        fit, s, _ = _mixed(y, C, g)
        cand["linear"] = (fit.bic, *_effect_at(fit, s, cmax, ndf), np.nan)

        grid = np.quantile(C[C > 0], np.linspace(0.05, 0.99, 40))
        for name, tf in (("log-linear", lambda c, p: np.log1p(c / p)),
                         ("Emax", lambda c, p: c / (c + p))):
            best = None
            for p in grid:
                f, ss, _ = _mixed(y, tf(C, p), g)
                if best is None or f.bic < best[0]:
                    best = (f.bic, *_effect_at(f, ss, tf(cmax, p), ndf), p)
            cand[name] = best

        bics = {k: v[0] for k, v in cand.items()}
        best_name = min(bics, key=bics.get)
        for name, (bic, eff, lo, hi, par) in cand.items():
            rows.append(dict(arm="in vivo", drug=drug, model=name, bic=bic,
                             delta_bic=bic - bics[best_name], is_best=(name == best_name),
                             effect_at_cmax=eff, ci_lo=lo, ci_hi=hi,
                             exceeds_10ms=bool(hi > 10), shape_param=par, n=len(on)))

    # the in-silico arm fits nothing: the simulation returns the curve
    od = pd.read_csv(RES / "ord_apd90.csv")
    o = od[(od.drug == "dofetilide") & (od.cell == "mid") & od.apd90.notna() & (od.v_rest < -60)]
    lin = np.polyfit(o.multiple, o.dAPD90, 1)
    resid = o.dAPD90 - np.polyval(lin, o.multiple)
    rows.append(dict(arm="in silico", drug="Dofetilide", model="no model is fitted",
                     bic=np.nan, delta_bic=np.nan, is_best=True, effect_at_cmax=np.nan,
                     ci_lo=np.nan, ci_hi=np.nan, exceeds_10ms=None,
                     shape_param=float(np.max(np.abs(resid))), n=len(o)))

    out = pd.DataFrame(rows)
    out.to_csv(OUT, index=False)
    with pd.option_context("display.width", 200):
        print(out.round(2).to_string(index=False))
    print("\nshape_param: EC50 (ng/mL) for Emax, p2 for log-linear;")
    print("for the in-silico row it is the largest residual (ms) a straight line would leave.")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
