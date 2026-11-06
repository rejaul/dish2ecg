
"""O'Hara-Rudy 2011 multichannel drug block -> APD90 and transmural repolarisation.

Model : myokit models repo, c/ohara-2011.mmt (O'Hara et al. 2011, PLoS Comput Biol).
Block : Hill,  scale = 1 - C^h / (IC50^h + C^h),  applied to the conductance BASE of
        each channel so the model's own cell-type factors still apply.
IC50s : CiPA optimal Hill fits (FDA/CiPA repo).
"""
import os, sys, numpy as np, pandas as pd, myokit
from multiprocessing import Pool

CH = {"hERG":"ikr.gKr", "IKs":"iks.gKs", "IK1":"ik1.gK1", "Ito":"ito.gto",
      "INa":"ina.gNa", "INaL":"inal.gNaL", "ICaL":"ical.PCa_base"}
CELLS = {"endo":0, "epi":1, "mid":2}
CL, PREPACE, DRUGPACE = 1000.0, 1000, 300
LOGT = np.arange(0.0, CL, 0.05)

def build(mode, scales=None):
    m, _, _ = myokit.load("ohara-2011.mmt")
    m.get("cell.mode").set_rhs(myokit.Number(mode))
    for ch, path in (scales or {}).items():
        pass
    if scales:
        for path, s in scales.items():
            v = m.get(path)
            v.set_rhs(myokit.Multiply(myokit.Number(float(s)), v.rhs()))
    p = myokit.pacing.blocktrain(period=CL, duration=0.5, offset=50.0)
    return myokit.Simulation(m, p)

def apd90(t, v):
    t = np.asarray(t, float); v = np.asarray(v, float)
    dv = np.gradient(v, t)
    i_up = int(np.argmax(dv))
    t_up = t[i_up]
    i_pk = i_up + int(np.argmax(v[i_up:]))
    v_pk = v[i_pk]
    v_rest = float(np.min(v[:i_up + 1]))
    v90 = v_pk - 0.9 * (v_pk - v_rest)
    idx = np.where(v[i_pk:] <= v90)[0]
    if len(idx) == 0:
        return np.nan, np.nan, v_pk, v_rest
    i90 = i_pk + int(idx[0])
    if i90 == 0:
        t90 = t[0]
    else:
        t90 = float(np.interp(v90, [v[i90], v[i90 - 1]], [t[i90], t[i90 - 1]]))
    return float(t_up), float(t90), float(v_pk), float(v_rest)

def steady(mode):
    s = build(mode)
    s.pre(PREPACE * CL)
    return s.state()

def run_one(args):
    cell, mode, state0, drug, mult, conc, scales = args
    try:
        s = build(mode, scales)
        s.set_state(state0)
        s.pre(DRUGPACE * CL)
        d = s.run(CL, log=["engine.time", "membrane.V"], log_times=LOGT)
        t_up, t90, v_pk, v_rest = apd90(d["engine.time"], d["membrane.V"])
        return dict(cell=cell, drug=drug, multiple=mult, conc_nM=conc,
                    apd90=t90 - t_up if t90 == t90 else np.nan,
                    t_rep90=t90, v_peak=v_pk, v_rest=v_rest, ok=True)
    except Exception as e:
        return dict(cell=cell, drug=drug, multiple=mult, conc_nM=conc,
                    apd90=np.nan, t_rep90=np.nan, v_peak=np.nan, v_rest=np.nan,
                    ok=False, err=f"{type(e).__name__}: {e}")

def scales_for(ic, drug, conc_nM):
    out = {}
    sub = ic[ic.drug == drug]
    for r in sub.itertuples():
        if r.channel not in CH or not np.isfinite(r.ic50_nM) or not np.isfinite(r.hill_coefficient):
            continue
        blk = conc_nM ** r.hill_coefficient / (r.ic50_nM ** r.hill_coefficient + conc_nM ** r.hill_coefficient)
        out[CH[r.channel]] = 1.0 - blk
    return out

if __name__ == "__main__":
    nproc = int(sys.argv[1]) if len(sys.argv) > 1 else 16
    ic = pd.read_csv("ic50.csv"); cond = pd.read_csv("conditions.csv")

    print("pre-pacing controls ...", flush=True)
    base = {c: steady(mo) for c, mo in CELLS.items()}
    ctl = [run_one((c, CELLS[c], base[c], "control", 0.0, 0.0, None)) for c in CELLS]
    C = pd.DataFrame(ctl)
    print(C[["cell", "apd90", "t_rep90", "v_peak", "v_rest"]].to_string(index=False), flush=True)

    # --- sanity gate: a 50% ICaL block MUST shorten APD90 ---
    t = run_one(("epi", 1, base["epi"], "_test", 0, 0, {"ical.PCa_base": 0.5}))
    c0 = float(C.loc[C.cell == "epi", "apd90"].iloc[0])
    print(f"GATE ICaL x0.5 -> APD90 {t['apd90']:.1f} vs control {c0:.1f}", flush=True)
    assert t["apd90"] < c0 - 5, "ICaL block did not shorten APD90 - conductance scaling is not reaching the model"

    jobs = []
    for r in cond.itertuples():
        sc = scales_for(ic, r.drug, r.conc_nM)
        for c, mo in CELLS.items():
            jobs.append((c, mo, base[c], r.drug, r.multiple, r.conc_nM, sc))
    print(f"running {len(jobs)} simulations on {nproc} procs", flush=True)
    with Pool(nproc) as pool:
        res = pool.map(run_one, jobs, chunksize=4)
    R = pd.DataFrame(res)
    R = R.merge(C[["cell", "apd90", "t_rep90"]].rename(
        columns={"apd90": "apd90_ctl", "t_rep90": "t_rep90_ctl"}), on="cell", how="left")
    R["dAPD90"] = R.apd90 - R.apd90_ctl
    R.to_csv("ord_apd90.csv", index=False)

    # transmural summation (UNCOUPLED cells: upper bound on dispersion)
    g = R.pivot_table(index=["drug", "multiple"], columns="cell", values="t_rep90")
    g["qt_surrogate"] = g.max(axis=1)
    g["transmural_dispersion"] = g.max(axis=1) - g[["endo", "epi", "mid"]].min(axis=1)
    c_qt = float(C.t_rep90.max()); c_disp = float(C.t_rep90.max() - C.t_rep90.min())
    g["dQT_surrogate"] = g.qt_surrogate - c_qt
    g["d_dispersion"] = g.transmural_dispersion - c_disp
    g.reset_index().to_csv("ord_transmural.csv", index=False)
    print(f"control QT surrogate {c_qt:.1f} ms, dispersion {c_disp:.1f} ms", flush=True)
    print("failures:", int((~R.ok).sum()), flush=True)
    print("wrote ord_apd90.csv, ord_transmural.csv", flush=True)
