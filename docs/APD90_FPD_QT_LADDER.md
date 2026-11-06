# APD90 -> FPD -> QT: what of this can actually be harmonised?

An audit of the attached conversation transcript against the data, and an
implementation of the parts that survive the audit.

The chain in question:

    Single-cell APD90  ->  syncytial field potential FPD  ->  transmural & spatial summation (QT)

Short answer: **arrow 1 cannot be estimated from the only public multi-site
dataset, because it is perfectly confounded. Arrow 2 is already built. Arrow 3
is not a weighted average and does not survive the simplification the document
proposes.** Details below.

---

## 1. Audit of the attached document

The document's architecture is the standard one and matches what this repo
already does: the same rate-correction algebra on both sides, a free-Cmax
common scale, multichannel block rather than hERG alone, and an in-silico
model to bridge. Those framings are right. The following specifics are not.

| Claim in the document | What the data say |
|---|---|
| Median beat "improves the signal-to-noise ratio by roughly a factor of sqrt(N)" | That is the **mean's** property. For Gaussian noise the median is *less* efficient by sqrt(pi/2) ~ 1.25. The median is chosen for outlier and ectopic rejection, which the document also states correctly - not for SNR. |
| Prefer "dataset-specific linear regression correction slopes (FPD = a + b*ISI)" over Fridericia | Tested on the CiPA MEA pilot (Millard 2018), the only public source that publishes beat period, FPD **and** FPDc. Fridericia reproduces the published FPDc to a median of **0.6 ms**, 17 of 18 datasets within 5 ms, implied exponent **0.3244** vs 1/3. The standard exponent is adequate here, and the field uses a **power law**, not the linear form the document gives. |
| hiPSC-CM spontaneous rate "40-70 bpm" | Schocken 2018 measured interbeat intervals of **1643-2102 ms = 29-37 bpm**. The document's range is too fast. |
| Hysteresis: EWMA and bi-exponential models, tau fitted on Holter | Correct in general, inapplicable here. The public clinical waveforms are **10-second strips**, 8-17% of one QT adaptation time constant, so the determining RR history precedes the recording. Hysteresis is **structurally unidentifiable** in this corpus, not merely unimplemented (METHODS 3.8). |
| Normalise in-vitro by free fraction in medium (fu,medium) | Directionally unsafe. Schocken 2018 found **8 of 18 drugs induced arrhythmias at LOWER concentrations with serum present** - opposite to the protein-binding expectation, and the paper itself calls it counterintuitive. A naive fu correction would move results the wrong way for nearly half the panel. |
| "Empirical mapping function: dQTc = m * dAPD90,weighted + c" | A transmural *weighted average* is not what sets QT. In an uncoupled model the QT surrogate is just the **M-cell** APD90 in 121 of 136 conditions, including every therapeutic match point (section 4). Averaging is the wrong operator. |
| Minimal Myokit script: `m['ikr.gkr'] = gkr_baseline * (1-block)` | Does not run. The variable is `ikr.gKr` (case-sensitive) and assignment goes through `m.get(...).set_rhs(...)`. It also scales `gKr` only, while in the O'Hara-Rudy model `ical.PCa` is **derived** from `PCa_base` by a cell-type piecewise, and `PK`/`PNa` are proportional to `PCa` - so an ICaL block must be applied to the base, not the leaf. The working version is in `scripts/ord_simulation.py`. |

Not mentioned anywhere in the document, and decisive in practice: the
confounding in section 2, and informative missingness (wells that arrest are
dropped, and the drop rate rises from 0.9% to 60.8% across the concentration
range).

---

## 2. Arrow 1 (APD90 -> FPD) is not identifiable in public data

Blinova 2018 pools both endpoints into a **single column**, separated only by
platform: "baseline and vehicle-controlled, rate-corrected action potential
duration at 90% repolarization, ddAPD90c, or field potential duration, ddFPDc,
depending on the EP platform used". Four MEA platforms (AXN, ECR, MCS, AMD)
report FPD; exactly one voltage-sensing optical platform, CellOPTIQ (CLY),
reports APD90. The paper asserts the two are "equivalent measures".

**The optical arm is one site.** CLY appears only at site 2 - 1,344 wells, 672
iCell2 and 672 Cor.4U. And site 2 "was the only test site that used the VSO
platform instead of the MEA platforms **and serum-free experimental medium**
instead of serum-containing medium". Endpoint, site and culture medium vary
together.

Leave-one-site-out, each site's 28-drug x 4-concentration profile regressed on
the pooled mean of the other nine (112 cells per site):

| | Pearson r | gain vs the others |
|---|---|---|
| 9 MEA sites | 0.481 - 0.944 (median 0.858) | 0.230 - 1.103 (median 0.919) |
| **Site 2 (optical, APD90)** | **0.710** | **1.436 +- 0.146** |

The optical arm reads **~44% larger** than the consensus - the highest gain of
all ten sites, z = +2.16 against the MEA sites, CI excluding 1.0. Paired
directly against the MEA mean on the 97 drug x concentration cells both arms
measured: **optical = 1.51 x MEA + 24.8 ms** (95% CI 1.22-1.80, r = 0.723,
p = 6e-17). Per drug, the optical arm is larger in **21 of 28**, median
+47.8 ms at the top shared concentration.

So Blinova's "equivalent measures" does not hold as a **scale** equivalence.
But the 1.44 **cannot be attributed to the endpoint**: Schocken 2018 shows
medium alone shifts baseline FPDc by 104 ms (575 +- 59 with serum vs 471 +- 44
without), and its arrhythmia-threshold result runs opposite to the
protein-binding prediction - so the confound cannot even be signed, let alone
subtracted. **Arrow 1 is unidentifiable in public data.** Estimating it needs
one site running both readouts on the same plates in the same medium.

---

## 3. The in-silico route, which has no site and no serum

O'Hara-Rudy 2011, multichannel Hill block from the CiPA optimal IC50 fits
applied to the conductance bases, 12 drugs x 11-12 concentrations x 3 cell
types = **408 simulations, 0 failures**. Controls reproduce the published
model: APD90 endo **265.5**, epi **228.0**, mid **329.9 ms**, resting
potential -88 mV. A gate asserts a 50% ICaL block shortens APD90 (228.0 ->
216.0) so the block demonstrably reaches the model.

At the same observed free-Cmax match points used by the clinical-vs-dish bridge:

| Drug | x free Cmax | in-silico dQT | dish ddFPDc | clinical ddQTcF |
|---|---|---|---|---|
| Dofetilide | 1.074 | +50.4 | +132.0 | +73.1 |
| Quinidine | 0.288 | +45.2 | +135.1 | +74.9 |
| Ranolazine | 0.956 | +37.0 | +38.8 | +11.2 |
| Verapamil | 0.500 | +12.2 | **-86.3** | +3.5 |

**Sign agreement with the clinic: in-silico 4/4, dish 3/4.** The model given
verapamil's *own* channel profile prolongs slightly, as the clinic does; the
dish shortens by 86 ms. Since the simulation uses the same multichannel block
the dish cells are exposed to, verapamil's in-vitro failure is better read as
an assay property - ICaL dominance in immature, IK1-poor cells - than as a
property of the drug's channel profile.

Magnitude is another matter. Neither surrogate is calibrated: the in-silico
value is 0.60-0.69x the clinical effect for the two strong hERG blockers and
3.3-3.5x for the two multichannel drugs. It is a sign-and-rank instrument, not
a millisecond predictor.

---

## 4. Arrow 3 is not a weighted average

With uncoupled cells the "QT surrogate" - the last repolarisation time across
the wall - is simply the **M-cell APD90**: the M cell repolarises last in **121
of 136 conditions**, including every therapeutic match point. (The 15
exceptions are all at **>= 3.78x free Cmax**, where the M cell loses a stable
single repolarisation and two fail outright - the M cell destabilises first,
which is itself the expected proarrhythmic signature.) The implied transmural
dispersion is **101.9 ms**, far above the ~30-40 ms a real coupled wall shows. Electrotonic coupling is the entire content of arrow 3: it is what
pulls the M cell back toward its neighbours. A weighted average of three
uncoupled cells cannot represent it, so the document's empirical mapping
function has no mechanism behind it.

**The coupled 1D cable was attempted and is blocked on this host.** A 165-cell
transmural wedge (60 endo / 45 M / 60 epi, dx = 0.01 cm) was built and the
model patched with the `diffusion_current` binding that a single-cell model
file lacks - without it myokit integrates every cell independently and
`set_conductance()` silently does nothing, which is a trap worth recording.
It then fails for a second reason: the available OpenCL runtime (pocl)
advertises `cl_khr_fp64` on the **device** but not on the **platform**, myokit
therefore refuses double precision, and O'Hara-Rudy is not numerically stable
in single precision under explicit Euler (NaN at dt = 0.002 ms). Completing
this rung needs a vendor OpenCL ICD with platform-level fp64, or the CiPA
pipeline's own 1D code. `scripts/cable_simulation.py` is left ready to run.

---

## 5. What this means for the bridge

The repo's position is unchanged and slightly better supported. The in-vitro
endpoint still must not be expressed as a QT interval in milliseconds: three
independent reasons now, not two - the dish-clinic scale gap, the
APD90-vs-FPD gain of ~1.44 that cannot be resolved, and the fact that the
transmural step is a coupling problem rather than an algebraic one.

What the simulation adds is a **confound-free sign check**. Where dish and
clinic disagree, running the drug's own channel profile through a human
ventricular model says whether the disagreement is in the drug or in the
assay. For verapamil it is the assay.
