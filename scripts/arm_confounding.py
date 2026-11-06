"""Is the APD90 arm of Blinova 2018 separable from site and medium?  (No.)

Blinova 2018 pools ddAPD90c and ddFPDc into ONE column, split only by platform:
four MEA platforms (AXN, ECR, MCS, AMD) report FPD, and one voltage-sensing
optical platform, CellOPTIQ (CLY), reports APD90. CLY is used by site 2 only,
which is also the only serum-free site - so endpoint, site and medium are
perfectly confounded.  This script quantifies the gain and shows it sits
outside the MEA site-to-site spread, which bounds but does not identify it.
"""
import argparse, pathlib, numpy as np, pandas as pd
from scipy import stats


def site_profiles(wells):
    m = wells.dropna(subset=["ddFPDc"])
    return m.pivot_table(index="site", columns=["drug", "conc"],
                         values="ddFPDc", aggfunc="mean"), m


def leave_one_site_out(prof, m):
    rows = []
    for s in prof.index:
        ref = prof.drop(index=s).mean(axis=0)
        ok = prof.loc[s].notna() & ref.notna()
        x, y = ref[ok].values, prof.loc[s][ok].values
        lr = stats.linregress(x, y)
        rows.append(dict(site=s, platform=m.loc[m.site == s, "Platform"].mode()[0],
                         arm="VSO (APD90)" if s == 2 else "MEA (FPD)",
                         n_cells=int(ok.sum()), pearson_r=stats.pearsonr(x, y)[0],
                         gain=lr.slope, gain_se=lr.stderr,
                         mean_signed_diff=float(np.mean(y - x))))
    return pd.DataFrame(rows).sort_values("pearson_r")


def paired_arms(m):
    m = m.assign(arm=np.where(m.Platform == "CLY", "VSO", "MEA"))
    t = m.pivot_table(index=["drug", "conc"], columns="arm",
                      values="ddFPDc", aggfunc="mean").dropna()
    lr = stats.linregress(t.MEA, t.VSO)
    return t, lr


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--wells", required=True)
    ap.add_argument("--out", default="results")
    a = ap.parse_args()
    out = pathlib.Path(a.out); out.mkdir(parents=True, exist_ok=True)
    wells = pd.read_csv(a.wells)
    prof, m = site_profiles(wells)
    L = leave_one_site_out(prof, m)
    L.to_csv(out / "loso_site_profiles.csv", index=False)
    mea, v = L[L.arm == "MEA (FPD)"], L[L.site == 2].iloc[0]
    print(f"MEA sites: r {mea.pearson_r.min():.3f}-{mea.pearson_r.max():.3f} "
          f"(median {mea.pearson_r.median():.3f}), gain {mea.gain.min():.3f}-"
          f"{mea.gain.max():.3f} (median {mea.gain.median():.3f})")
    print(f"VSO site 2: r {v.pearson_r:.3f}, gain {v.gain:.3f} +- {v.gain_se:.3f}  "
          f"z(gain) vs MEA = {(v.gain-mea.gain.mean())/mea.gain.std():+.2f}")
    t, lr = paired_arms(m)
    t.reset_index().to_csv(out / "arm_bycell.csv", index=False)
    # per-drug view at the highest concentration both arms measured
    per = []
    for d, g in t.groupby(level=0):
        g = g.reset_index()
        top = g.iloc[g["conc"].idxmax()]
        per.append(dict(drug=d, conc_idx=int(top["conc"]), MEA=top["MEA"], VSO=top["VSO"],
                        diff=top["VSO"] - top["MEA"],
                        ratio=top["VSO"] / top["MEA"] if abs(top["MEA"]) > 5 else np.nan))
    pd.DataFrame(per).sort_values("diff").to_csv(out / "arm_perdrug.csv", index=False)
    print(f"optical = {lr.slope:.3f} x MEA + {lr.intercept:.1f}  "
          f"(95% CI {lr.slope-1.96*lr.stderr:.3f}-{lr.slope+1.96*lr.stderr:.3f}, "
          f"r={lr.rvalue:.3f}, n={len(t)})")
    print("NOTE: the gain is NOT attributable to the endpoint - site 2 is also the "
          "only serum-free site, and Schocken 2018 shows medium alone moves baseline "
          "FPDc by 104 ms in an unpredictable direction.")


if __name__ == "__main__":
    main()
