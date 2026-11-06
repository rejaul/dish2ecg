# Does the simulation sort drugs by TdP risk better than the dish?

> **SUPERSEDED IN PART — see `docs/RISK_CLASSIFICATION_28.md`.** This document
> covers the 12 CiPA *training* drugs only. The full 28-drug test, with 16 drugs
> held out, has since been run. Its result reverses this document's conclusion
> about the dish: the claim that the dish "detects danger but does not grade
> safety" is an artifact of the 12-drug subset and is withdrawn.

The dish exists to assign drugs to high / intermediate / low torsades risk.
This compares the CiPA in-silico metric (qNet) against the hiPSC-CM endpoints on
exactly the same 12 drugs.

**The headline caveat first, because it governs everything below: those 12 drugs
are the CiPA *training* set.** The ORd-CiPA v1.0 model and the qNet metric were
optimised on precisely these compounds (Dutta 2017). qNet's performance here is
therefore **in-sample**, and this is *not* an independent test. What it is: a
faithful reproduction of the published training result, which validates the
implementation, plus an out-of-sample look at the dish on the same drugs. **No
claim that the simulation beats the dish is supported by this evidence.**

---

## 1. Implementation, and how it was checked

ORd-CiPA v1.0 (Dutta 2017), endocardial cell, 0.5 Hz, doses 1-4x therapeutic -
the CiPA protocol. IKr block uses the model's **dynamic Markov binding**
(drug-specific Kmax, Ku, n, halfmax, Vhalf from the CiPA hERG fits), not
conductance scaling; the other six channels use Hill block from the CiPA optimal
IC50 fits. 48 drug-dose simulations, 0 failures.

Three implementation details cost real time and are easy to get silently wrong:

| | |
|---|---|
| **qNet's current set** | `newordherg_qNet.c` line 682 is `ydot[49] = INaL+ICaL+Ito+IKr+IKs+IK1`. **INa is excluded**, and their `ICaL` is the **calcium carrier alone** (ICaNa and ICaK are separate variables). Including INa put control qNet at 0.0044 instead of 0.0785 uC/uF - an 18-fold error that still looks like a plausible number. |
| **Kt** | A fixed, non-drug-specific rate (`newordherg_pars.txt`: `Kt 3.5e-05`) that the myokit model file defaults to **0**, switching off the Cbound trapping transitions. With Kt = 0, dofetilide at 4x prolongs APD90 by only 28 ms; with Kt set, 130 ms. |
| **Drug concentration units** | The mmt writes `r1` with `(D / 1 [mM])`, so D is used dimensionlessly, while CiPA calibrates halfmax in nM^n (`AP_simulation.R`: "cmax ... should be in nanomolar"). `ikr.D` must carry the **nanomolar number**. A gate in the script fails loudly if this is wrong - setting D in mM produces no block at all. |

Checks that passed: control qNet **0.0785 uC/uF** and control APD90 **286.7 ms**,
both consistent with published ORd-CiPA values; dofetilide at 4x therapeutic
prolongs APD90 to **416.8 ms** and drops qNet to 0.0523.

## 2. Result

qNet puts the three risk categories in **completely non-overlapping bands**:

| CiPA risk | qNet range (uC/uF) |
|---|---|
| low | 0.0831 - 0.1009 |
| intermediate | 0.0687 - 0.0747 |
| high | 0.0390 - 0.0648 |

which reproduces the published training-set separation - the point of running it.

Against the dish endpoints on the same 12 drugs:

| predictor | Spearman rho with risk | AUC high vs rest | AUC low vs rest |
|---|---|---|---|
| qNet (in-silico) | **0.946** [0.85, 0.96] | 1.000 | 1.000 |
| arrhythmia-event fraction (dish) | 0.650 [0.06, 0.90] | **1.000** | 0.688 [0.36, 1.00] |
| max ddFPDc (dish) | 0.680 [0.27, 0.89] | 0.844 [0.56, 1.00] | 0.875 [0.60, 1.00] |
| ddFPDc at Cmax (dish) | 0.680 [0.13, 0.93] | 0.938 [0.72, 1.00] | 0.781 [0.40, 1.00] |

95% CIs are 4,000-sample bootstraps over the 12 drugs. Where separation is
complete the bootstrap CI collapses to [1.000, 1.000] and is **uninformative** -
it reflects the resampling, not the uncertainty.

The informative finding is not the qNet column, which is in-sample. It is the
**shape of the dish's competence**: its arrhythmia-event fraction identifies all
four high-risk drugs perfectly (AUC 1.00) - one of the three predictors Blinova 2018 found
most useful, alongside prolongation at any concentration and at clinical Cmax - while separating low from the rest at only AUC 0.69,
with a CI running down to 0.36. The dish's own prolongation endpoints are weaker
at both ends than its event counts. **The dish detects danger; it does not grade
safety.**

## 3. What would make this a real test

The CiPA *validation* set is 16 further drugs, and all 16 are in the Blinova
28-drug panel - so the dish side is already in hand. What is missing is their
ion-channel and hERG-binding parameters: the FDA/CiPA repository publishes
Hill and hERG fits for the 12 training drugs only. Those parameters are in
Li et al. 2019's supplement. With them, the same pipeline runs unchanged and
the comparison becomes out-of-sample on both arms.
