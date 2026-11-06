
"""Out-of-sample CiPA qNet: all 28 drugs, hERG uncertainty by bootstrap ensemble.

Model   : ORd-CiPA v1.0 (Dutta 2017), endocardial, 0.5 Hz, doses 1-4x therapeutic.
IKr     : dynamic Markov binding. The FDA repo publishes an OPTIMAL point fit only
          for the 12 training drugs; for all 28 it publishes 2000 nonparametric
          BOOTSTRAP parameter sets. Marginal medians of those samples are not a
          safe proxy (bepridil's Kmax differs 5-fold from the point fit), because
          Kmax and halfmax are strongly co-dependent. So qNet is computed for each
          bootstrap sample and summarised in OUTPUT space: median qNet is the point
          estimate, percentiles give the interval. The 12 training drugs, where the
          published point fit exists, validate that choice.
Others  : Hill block from newCiPA.csv. IC50 <= 0 or hill <= 0 means the channel was
          not measured -> NO block (naively applying the Hill formula to IC50 = 0
          would give 100% block).
qNet    : integral over one beat of INaL+ICaL+Ito+IKr+IKs+IK1, in uC/uF
          (FDA newordherg_qNet.c line 682; INa excluded, ICaL = calcium carrier).
"""
import sys, numpy as np, pandas as pd, myokit
from multiprocessing import Pool

CURRENTS = ["inal.INaL", "ical.ICaLCa", "ito.Ito", "ikr.IKr", "iks.IKs", "ik1.IK1"]
CL, PREPACE, KT = 2000.0, 1000, 3.5e-5
LOGT = np.arange(0.0, CL, 0.05)


def build(herg=None, hill=None):
    m, _, _ = myokit.load("ohara-cipa-v1-2017.mmt")
    m.get("cell.mode").set_rhs(myokit.Number(0))
    for k, v in (herg or {}).items():
        m.get(f"ikr.{k}").set_rhs(myokit.Number(float(v)))
    for path, s in (hill or {}).items():
        v = m.get(path); v.set_rhs(myokit.Multiply(myokit.Number(float(s)), v.rhs()))
    p = myokit.pacing.blocktrain(period=CL, duration=0.5, offset=50.0)
    return myokit.Simulation(m, p)


def qnet_apd(sim, state=None):
    if state is not None:
        sim.set_state(state)
    sim.pre(PREPACE * CL)
    d = sim.run(CL, log=["engine.time", "membrane.V"] + CURRENTS, log_times=LOGT)
    t = np.asarray(d["engine.time"], float)
    q = float(np.trapezoid(sum(np.asarray(d[c], float) for c in CURRENTS), t) / 1000.0)
    v = np.asarray(d["membrane.V"], float)
    dv = np.gradient(v, t); i_up = int(np.argmax(dv))
    i_pk = i_up + int(np.argmax(v[i_up:]))
    v90 = v[i_pk] - 0.9 * (v[i_pk] - float(np.min(v[:i_up + 1])))
    idx = np.where(v[i_pk:] <= v90)[0]
    a = float(t[i_pk + idx[0]] - t[i_up]) if len(idx) else np.nan
    return q, a, sim.state()


def hill_scales(H, drug, conc_nM):
    out = {}
    for r in H[H.drug == drug].itertuples():
        if not (np.isfinite(r.ic50_nM) and np.isfinite(r.hill)) or r.ic50_nM <= 0 or r.hill <= 0:
            continue                                   # channel not measured -> no block
        blk = conc_nM ** r.hill / (r.ic50_nM ** r.hill + conc_nM ** r.hill)
        out[r.target] = 1.0 - blk
    return out


def one(a):
    drug, samp, dose, conc, herg, hill, st0 = a
    try:
        q, ap, _ = qnet_apd(build(herg, hill), st0)
        return dict(drug=drug, sample=samp, dose=dose, conc_nM=conc, qnet=q, apd90=ap, ok=True)
    except Exception as e:
        return dict(drug=drug, sample=samp, dose=dose, conc_nM=conc, qnet=np.nan,
                    apd90=np.nan, ok=False, err=f"{type(e).__name__}: {e}")


if __name__ == "__main__":
    nproc = int(sys.argv[1]) if len(sys.argv) > 1 else 24
    V = pd.read_csv("newcipa28.csv"); H = pd.read_csv("hill28.csv"); B = pd.read_csv("herg_boot28.csv")
    q0, a0, st0 = qnet_apd(build())
    print(f"CONTROL qNet {q0:.5f} uC/uF, APD90 {a0:.1f} ms", flush=True)
    assert 0.07 < q0 < 0.09, f"control qNet {q0} does not match the established control"

    ther = dict(zip(V.dl, V.therapeutic))
    jobs = []
    for r in B.itertuples():
        h = {"Kmax": r.Kmax, "Ku": r.Ku, "n": r.n, "halfmax": r.halfmax, "Vhalf": r.Vhalf, "Kt": KT}
        for dose in (1, 2, 3, 4):
            c = dose * ther[r.drug]
            jobs.append((r.drug, r.sample, dose, c, {**h, "D": c}, hill_scales(H, r.drug, c), st0))
    print(f"running {len(jobs)} simulations on {nproc} procs", flush=True)
    with Pool(nproc) as pool:
        res = pool.map(one, jobs, chunksize=2)
    R = pd.DataFrame(res); R["qnet_control"] = q0
    R.to_csv("qnet28_raw.csv", index=False)
    print("failures:", int((~R.ok).sum()), flush=True)

    per = R.groupby(["drug", "sample"]).qnet.mean().reset_index(name="qnet_1to4")
    S = (per.groupby("drug").qnet_1to4
           .agg(qnet_median="median", qnet_lo=lambda s: s.quantile(0.025),
                qnet_hi=lambda s: s.quantile(0.975), n_samples="size")
           .reset_index().merge(V[["dl", "drug", "CiPA", "split", "therapeutic"]]
                                .rename(columns={"drug": "drug_name", "dl": "drug"}), on="drug"))
    S.to_csv("qnet28_by_drug.csv", index=False)
    print(S.sort_values("qnet_median").to_string(index=False), flush=True)
