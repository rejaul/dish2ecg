
"""1D transmural cable (O'Hara-Rudy 2011) with pseudo-ECG.

Wedge  : 60 endo | 45 M | 60 epi, dx = 0.01 cm  (1.65 cm wall), stimulated at the
         endocardial end. Cell type is imposed as a FIELD on cell.mode, so every
         mode-dependent parameter in the model follows automatically.
Drug   : Hill block applied uniformly to the conductance bases (same as single cell).
PseudoECG: Gima & Rudy 2002, electrode 2.0 cm beyond the epicardial end;
         phi(t) ~ sum_i (-dV/dx)_i * (1/r_i^2) * dx   (timing only, scale dropped).
"""
import os, sys, json, numpy as np, pandas as pd, myokit

CH = {"hERG":"ikr.gKr","IKs":"iks.gKs","IK1":"ik1.gK1","Ito":"ito.gto",
      "INa":"ina.gNa","INaL":"inal.gNaL","ICaL":"ical.PCa_base"}
N_ENDO,N_M,N_EPI = 60,45,60
NCELL = N_ENDO+N_M+N_EPI
DX = 0.01                      # cm
GX = 11.0                      # mS/uF  -> D = GX*dx^2 = 1.1e-3 cm^2/ms
CL = 1000.0
ELECTRODE_X = NCELL*DX + 2.0   # cm
MODES = np.array([0]*N_ENDO + [2]*N_M + [1]*N_EPI, dtype=float)   # endo=0, epi=1, mid=2

def _select_device():
    """Pick an OpenCL device that really supports double precision."""
    import pyopencl as cl
    for p_ in cl.get_platforms():
        for d_ in p_.get_devices():
            if "cl_khr_fp64" in d_.extensions:
                myokit.OpenCL.save_selection(p_.name, d_.name)
                print("OpenCL device:", p_.name, "/", d_.name, flush=True)
                return
    raise RuntimeError("no OpenCL device with cl_khr_fp64")

_select_device()

def scaled_model(scales, cable=False):
    m,_,_ = myokit.load("ohara-2011.mmt")
    for path, s in (scales or {}).items():
        v = m.get(path); v.set_rhs(myokit.Multiply(myokit.Number(float(s)), v.rhs()))
    if cable:
        # A single-cell model file has no diffusion_current binding, so every cell
        # would integrate independently and set_conductance() would do nothing.
        mem = m.get("membrane")
        idiff = mem.add_variable("i_diff")
        idiff.set_rhs(myokit.Number(0))
        idiff.set_unit(myokit.parse_unit("A/F"))
        idiff.set_binding("diffusion_current")
        v = m.get("membrane.V")
        v.set_rhs(myokit.Minus(v.rhs(), myokit.Name(idiff)))   # dot(V) = -(i_ion+i_stim+i_diff)
    return m

def single_cell_steady(mode, scales, beats=1000):
    m = scaled_model(scales)
    m.get("cell.mode").set_rhs(myokit.Number(mode))
    p = myokit.pacing.blocktrain(period=CL, duration=0.5, offset=50.0)
    s = myokit.Simulation(m, p); s.pre(beats*CL)
    return s.state()

def t90_of(t, v):
    t=np.asarray(t,float); v=np.asarray(v,float)
    dv=np.gradient(v,t); i_up=int(np.argmax(dv))
    i_pk=i_up+int(np.argmax(v[i_up:])); v_pk=v[i_pk]
    v_rest=float(np.min(v[:i_up+1])); v90=v_pk-0.9*(v_pk-v_rest)
    idx=np.where(v[i_pk:]<=v90)[0]
    if len(idx)==0: return np.nan, np.nan
    i90=i_pk+int(idx[0])
    t90 = float(np.interp(v90,[v[i90],v[i90-1]],[t[i90],t[i90-1]])) if i90>0 else t[i90]
    return float(t[i_up]), t90

def pseudo_ecg(V, dx=DX, xe=ELECTRODE_X):
    # V: (ntime, ncell)
    x = (np.arange(V.shape[1])+0.5)*dx
    r = xe - x
    dVdx = np.gradient(V, dx, axis=1)
    return -(dVdx * (1.0/r**2)).sum(axis=1)*dx

def tangent_tend(t, phi):
    """T-end by max-downslope tangent to baseline, on the T wave."""
    i_t = int(np.argmax(np.abs(phi[len(phi)//4:]))) + len(phi)//4   # T peak (largest |phi| late)
    seg = slice(i_t, len(phi))
    d = np.gradient(phi[seg], t[seg])
    i_s = i_t + int(np.argmax(np.abs(d)) if phi[i_t] > 0 else np.argmax(np.abs(d)))
    slope = np.gradient(phi, t)[i_s]
    if slope == 0: return np.nan
    return float(t[i_s] - phi[i_s]/slope)

def run_cable(scales, label, beats=3):
    states = {mo: single_cell_steady(mo, scales) for mo in (0,1,2)}
    m = scaled_model(scales, cable=True)
    p = myokit.pacing.blocktrain(period=CL, duration=0.5, offset=50.0)
    # pocl advertises cl_khr_fp64 on the DEVICE but not the PLATFORM, so myokit
    # falls back to single precision; an explicit small Euler step is then required.
    s = myokit.SimulationOpenCL(m, p, ncells=NCELL, diffusion=True,
                                precision=myokit.SINGLE_PRECISION)
    s.set_step_size(0.002)
    s.set_conductance(GX)
    s.set_paced_cells(nx=5)                      # stimulate the endocardial end
    s.set_field(m.get("cell.mode"), MODES.tolist())
    full = []
    for i, mo in enumerate(MODES.astype(int)):
        full.extend(states[mo])
    s.set_state(full)
    s.run(beats*CL - CL, log=myokit.LOG_NONE)
    d = s.run(CL, log=["engine.time","membrane.V"], log_interval=0.1)
    t = np.asarray(d["engine.time"], float)
    V = np.column_stack([d[f"{i}.membrane.V"] for i in range(NCELL)])
    n_act = int((V.max(axis=0) > 0).sum())
    if n_act < 0.9*NCELL:
        raise RuntimeError(f"conduction failure in '{label}': only {n_act}/{NCELL} cells depolarised")
    tup = np.full(NCELL, np.nan); t90 = np.full(NCELL, np.nan)
    for i in range(NCELL):
        tup[i], t90[i] = t90_of(t, V[:,i])
    act = np.nanmin(tup)
    cv = (NCELL-10)*DX / ((np.nanmax(tup)-np.nanmin(tup))/1000.0) if np.nanmax(tup)>np.nanmin(tup) else np.nan
    phi = pseudo_ecg(V)
    t_end = tangent_tend(t, phi)
    return dict(label=label, act_time=act,
                apd90_endo=float(np.nanmean(t90[:N_ENDO]-tup[:N_ENDO])),
                apd90_mid =float(np.nanmean(t90[N_ENDO:N_ENDO+N_M]-tup[N_ENDO:N_ENDO+N_M])),
                apd90_epi =float(np.nanmean(t90[N_ENDO+N_M:]-tup[N_ENDO+N_M:])),
                disp_rep=float(np.nanmax(t90)-np.nanmin(t90)),
                last_rep=float(np.nanmax(t90)), cv_cm_s=float(cv),
                qt_pseudo=float(t_end-act) if t_end==t_end else np.nan), t, phi, V, tup, t90

def scales_for(ic, drug, conc_nM):
    out={}
    for r in ic[ic.drug==drug].itertuples():
        if r.channel not in CH or not np.isfinite(r.ic50_nM) or not np.isfinite(r.hill_coefficient): continue
        blk = conc_nM**r.hill_coefficient/(r.ic50_nM**r.hill_coefficient + conc_nM**r.hill_coefficient)
        out[CH[r.channel]] = 1.0-blk
    return out

if __name__ == "__main__":
    ic = pd.read_csv("ic50.csv"); cond = pd.read_csv("cable_conditions.csv")
    rows=[]; traces={}
    res, t, phi, V, tup, t90 = run_cable(None, "control")
    print("CONTROL", {k: (round(v,1) if isinstance(v,float) else v) for k,v in res.items()}, flush=True)
    assert 20 < res["cv_cm_s"] < 120, f"conduction velocity {res['cv_cm_s']:.0f} cm/s is non-physiological"
    rows.append(dict(drug="control", multiple=0.0, conc_nM=0.0, **res))
    traces["control"] = dict(t=t.tolist(), phi=phi.tolist())
    np.save("cable_V_control.npy", V.astype(np.float32))
    for r_ in cond.itertuples():
        sc = scales_for(ic, r_.drug, r_.conc_nM)
        res, t, phi, V, tup, t90 = run_cable(sc, f"{r_.drug}@{r_.multiple}")
        print(res["label"], {k: (round(v,1) if isinstance(v,float) else v) for k,v in res.items()}, flush=True)
        rows.append(dict(drug=r_.drug, multiple=r_.multiple, conc_nM=r_.conc_nM, **res))
        traces[res["label"]] = dict(t=t.tolist(), phi=phi.tolist())
    R = pd.DataFrame(rows)
    c0 = R[R.drug=="control"].iloc[0]
    R["dQT_pseudo"]  = R.qt_pseudo - c0.qt_pseudo
    R["d_disp"]      = R.disp_rep  - c0.disp_rep
    R["dLastRep"]    = R.last_rep  - c0.last_rep
    R.to_csv("cable_results.csv", index=False)
    json.dump(traces, open("cable_traces.json","w"))
    print("wrote cable_results.csv", flush=True)
