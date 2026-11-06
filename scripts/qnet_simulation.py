
"""CiPA qNet: ORd-CiPA v1.0 (Dutta 2017) with dynamic hERG binding + Hill block.

Model  : myokit models repo, c/ohara-cipa-v1-2017.mmt (Dutta et al. 2017, built
         from FDA/CiPA AP_simulation/models/newordherg_qNet.c).
IKr    : dynamic drug binding via the Markov model's own Kmax/Ku/n/halfmax/Vhalf,
         NOT conductance scaling.  Parameters are the CiPA optimal hERG fits.
Others : Hill conductance block from the CiPA optimal IC50 fits (hERG excluded).
Protocol: endocardial cell, 0.5 Hz (CL = 2000 ms), doses 1-4x therapeutic, as CiPA.
qNet   : integral over one beat of INet = INaL+ICaL+Ito+IKr+IKs+IK1,
         in uC/uF  ( = integral of [A/F] dt[ms] / 1000 ).

UNITS NOTE: the mmt writes r1 with (D / 1 [mM]), so D's value is used
dimensionlessly.  CiPA calibrates halfmax in nM^n (AP_simulation.R: "cmax ...
should be in nanomolar"), so ikr.D is set to the NUMERIC nanomolar value.
A gate below fails loudly if that convention is wrong.
"""
import sys, numpy as np, pandas as pd, myokit
from multiprocessing import Pool

HILL = {"ICaL":"ical.PCa_base", "INa":"ina.gNa", "INaL":"inal.gNaL",
        "Ito":"ito.gto", "IKs":"iks.gKs", "IK1":"ik1.gK1"}          # hERG excluded
# FDA/CiPA AP_simulation/models/newordherg_qNet.c line 682:
#   ydot[49] = INaL + ICaL + Ito + IKr + IKs + IK1
# INa is NOT in qNet, and their ICaL is the calcium carrier alone (ICaNa and
# ICaK are separate variables there) -> myokit's ical.ICaLCa.
CURRENTS = ["inal.INaL", "ical.ICaLCa", "ito.Ito", "ikr.IKr", "iks.IKs", "ik1.IK1"]
CL, PREPACE = 2000.0, 1000
# Kt is a FIXED, non-drug-specific rate in CiPA (newordherg_pars.txt: "Kt 3.5e-05").
# The mmt defaults it to 0, which switches off the Cbound trapping transitions
# entirely and badly under-states block for trapped drugs like dofetilide.
KT = 3.5e-5


def build(herg=None, hill=None):
    m, _, _ = myokit.load("ohara-cipa-v1-2017.mmt")
    m.get("cell.mode").set_rhs(myokit.Number(0))            # endocardial, as CiPA
    if herg:
        for k, v in herg.items():
            m.get(f"ikr.{k}").set_rhs(myokit.Number(float(v)))
    for path, s in (hill or {}).items():
        v = m.get(path); v.set_rhs(myokit.Multiply(myokit.Number(float(s)), v.rhs()))
    p = myokit.pacing.blocktrain(period=CL, duration=0.5, offset=50.0)
    return myokit.Simulation(m, p)


def qnet_and_apd(sim, state=None):
    if state is not None:
        sim.set_state(state)
    sim.pre(PREPACE * CL)
    d = sim.run(CL, log=["engine.time", "membrane.V"] + CURRENTS,
                log_times=np.arange(0.0, CL, 0.05))
    t = np.asarray(d["engine.time"], float)
    inet = sum(np.asarray(d[c], float) for c in CURRENTS)
    qnet = float(np.trapezoid(inet, t) / 1000.0)            # uC/uF
    v = np.asarray(d["membrane.V"], float)
    dv = np.gradient(v, t); i_up = int(np.argmax(dv))
    i_pk = i_up + int(np.argmax(v[i_up:]))
    v90 = v[i_pk] - 0.9 * (v[i_pk] - float(np.min(v[:i_up + 1])))
    idx = np.where(v[i_pk:] <= v90)[0]
    apd = float(t[i_pk + idx[0]] - t[i_up]) if len(idx) else np.nan
    return qnet, apd, sim.state()


def hill_scales(ic, drug, conc_nM):
    out = {}
    for r in ic[ic.drug == drug].itertuples():
        if r.channel not in HILL or not np.isfinite(r.ic50_nM) or not np.isfinite(r.hill_coefficient):
            continue
        blk = conc_nM ** r.hill_coefficient / (r.ic50_nM ** r.hill_coefficient + conc_nM ** r.hill_coefficient)
        out[HILL[r.channel]] = 1.0 - blk
    return out


def one(args):
    drug, dose, conc, herg, hill, state0 = args
    try:
        q, a, _ = qnet_and_apd(build(herg, hill), state0)
        return dict(drug=drug, dose=dose, conc_nM=conc, qnet=q, apd90=a, ok=True)
    except Exception as e:
        return dict(drug=drug, dose=dose, conc_nM=conc, qnet=np.nan, apd90=np.nan,
                    ok=False, err=f"{type(e).__name__}: {e}")


if __name__ == "__main__":
    nproc = int(sys.argv[1]) if len(sys.argv) > 1 else 16
    ic = pd.read_csv("ic50.csv")
    dr = pd.read_csv("cipa_drugs.csv")
    hb = pd.read_csv("herg_binding.csv").set_index("drug")

    q0, a0, st0 = qnet_and_apd(build())
    print(f"CONTROL qNet {q0:.5f} uC/uF, APD90 {a0:.1f} ms", flush=True)
    assert 0.01 < q0 < 0.25, f"control qNet {q0} outside the plausible CiPA range"

    # --- units gate: dofetilide 4x therapeutic (8 nM) must block hERG hard.
    hd = {**hb.loc["dofetilide"].to_dict(), "Kt": KT}
    q_nm, a_nm, _ = qnet_and_apd(build({**hd, "D": 8.0},
                                       hill_scales(ic, "dofetilide", 8.0)), st0)
    q_mm, a_mm, _ = qnet_and_apd(build({**hd, "D": 8.0e-6},
                                       hill_scales(ic, "dofetilide", 8.0)), st0)
    print(f"GATE dofetilide 4x: D as nM -> APD90 {a_nm:.1f} (qNet {q_nm:.5f}); "
          f"D as mM -> APD90 {a_mm:.1f} (qNet {q_mm:.5f}); control {a0:.1f}", flush=True)
    assert a_nm - a0 > 20, "nanomolar convention produced no hERG block - units are wrong"
    assert abs(a_mm - a0) < 5, "millimolar convention also acted - gate is not discriminating"

    jobs = []
    for r in dr.itertuples():
        h = {**hb.loc[r.drug].to_dict(), "Kt": KT}
        for dose in (1, 2, 3, 4):
            c = dose * r.therapeutic
            jobs.append((r.drug, dose, c, {**h, "D": c}, hill_scales(ic, r.drug, c), st0))
    print(f"running {len(jobs)} drug-dose simulations on {nproc} procs", flush=True)
    with Pool(nproc) as pool:
        res = pool.map(one, jobs, chunksize=1)
    R = pd.DataFrame(res)
    R["qnet_control"] = q0; R["apd90_control"] = a0
    R["dqnet"] = R.qnet - q0; R["dapd90"] = R.apd90 - a0
    R.to_csv("qnet_raw.csv", index=False)
    S = (R.groupby("drug").agg(qnet_mean_1to4=("qnet", "mean"), n_ok=("ok", "sum"))
           .reset_index().merge(dr, on="drug"))
    S.to_csv("qnet_by_drug.csv", index=False)
    print(f"failures: {int((~R.ok).sum())}", flush=True)
    print(S.sort_values("qnet_mean_1to4").to_string(index=False), flush=True)
