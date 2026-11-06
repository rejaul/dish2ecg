"""Does the published random-SLOPE specification close the ECGDMMLD dofetilide gap?

Johannesen L et al., Clin Pharmacol Ther 2015;99:214-223, doi:10.1002/cpt.205,
Statistical methods: "a linear mixed-effects model was used to explore the
concentration dependency of change in QTc, J-Tpeakc, and Tpeak-Tend with
concentration as a fixed effect and **subject as a random effect on intercept
and concentration**."

dish2ecg defaults to a random intercept only. This script refits both public
crossover studies both ways and reports the difference.
"""
from __future__ import annotations
import argparse, sys, warnings
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parent))
import clinical_reference as cr
from dish2ecg import fit_concentration_response
from dish2ecg.io import load_physionet_clinical

PUBLISHED = {("ecgdmmld", "Dofetilide", "QTcF"): 45.2,
             ("ecgdmmld", "Dofetilide", "JTpeakc"): 29.1,
             ("ecgdmmld", "Mexiletine + Dofetilide", "QTcF"): 17.3,
             ("ecgdmmld", "Mexiletine + Dofetilide", "JTpeakc"): 1.1,
             ("ecgrdvq", "Dofetilide", "QTcF"): 73.6,
             ("ecgrdvq", "Quinidine Sulph", "QTcF"): 78.9,
             ("ecgrdvq", "Ranolazine", "QTcF"): 12.0,
             ("ecgrdvq", "Dofetilide", "JTpeakc"): 39.1,
             ("ecgrdvq", "Quinidine Sulph", "JTpeakc"): 26.1}

def main(out: Path) -> None:
    warnings.filterwarnings("ignore")
    out.mkdir(parents=True, exist_ok=True)
    rows = []
    for study in ("ecgrdvq", "ecgdmmld"):
        raw = load_physionet_clinical(study)
        d, _ = cr.prepare(raw, study)
        g = cr.collapse_replicates(d).reset_index(drop=True)
        eps, dd, _, _ = cr.corrected_endpoints(g)
        for arm in sorted(g["arm"].unique()):
            if arm == "Placebo":
                continue
            sel = ((g["arm"] == arm) & g["conc_ngml"].notna() & (g["TPT"] > 0)).to_numpy()
            if sel.sum() < 10:
                continue
            cmax, _n = cr.gmean_cmax(g[sel])
            for nm in ("QTcF", "JTpeakc"):
                r = {"study": study, "arm": arm, "endpoint": nm,
                     "gmean_cmax": round(cmax, 3)}
                for lbl, rs in (("intercept_only", False), ("random_slope", True)):
                    try:
                        f = fit_concentration_response(
                            g[sel], dd[nm].subset(sel), concentration_col="conc_ngml",
                            group_col="RANDID", at_concentration=cmax, random_slope=rs)
                        r[lbl] = round(f.effect_at_concentration, 1)
                        r[lbl + "_lo"] = round(f.effect_ci90[0], 1)
                        r[lbl + "_hi"] = round(f.effect_ci90[1], 1)
                        r[lbl + "_boundary"] = bool(f.boundary_fallback)
                    except Exception as e:
                        r[lbl] = np.nan; r[lbl + "_err"] = repr(e)[:60]
                pub = PUBLISHED.get((study, arm, nm))
                r["published"] = pub
                if pub is not None:
                    r["gap_intercept_only"] = round(r["intercept_only"] - pub, 1)
                    r["gap_random_slope"] = round(r["random_slope"] - pub, 1)
                rows.append(r)
    T = pd.DataFrame(rows)
    T.to_csv(out / "random_slope_comparison.csv", index=False)
    cols = ["study", "arm", "endpoint", "intercept_only", "random_slope",
            "published", "gap_intercept_only", "gap_random_slope"]
    print(T[cols].to_string(index=False))
    v = T.dropna(subset=["published"])
    print("\n  mean |gap| intercept-only : %.2f ms" % v["gap_intercept_only"].abs().mean())
    print("  mean |gap| random-slope   : %.2f ms" % v["gap_random_slope"].abs().mean())

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--out", type=Path, default=Path("results"))
    main(ap.parse_args().out)
