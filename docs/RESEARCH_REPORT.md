# dish2ecg — research report

> **FROZEN AT 0.1.1 (86 tests).** This is the early technical reference. For current results read [`PROJECT_REPORT.md`](PROJECT_REPORT.md). Statements below that have since changed: (a) the missing fraction unbound for bepridil, chlorpromazine, cisapride and terfenadine was found (see `DRUG_REFERENCE_SOURCES.md`); (b) the T-end delineation rule was tested at beat level and does **not** explain the 6-7 ms ECGRDVQ-ECGCIPA shift (`DELINEATION_AND_SHIFT.md`); (c) the risk comparison was extended to all 28 drugs, with 16 held out (`RISK_CLASSIFICATION_28.md`), and the claim that the dish "does not grade safety" is withdrawn; (d) the waveforms, ORd simulations and optical-versus-electrode analysis post-date this version.


**Corrected repolarisation endpoints from hiPSC-CM recordings and from clinical ECG, computed by one cascade and compared on one concentration scale.**

Version 0.1.1 · 4 November 2026 · 86 tests · every number below is reproduced by a script in this repository. Citations are written `[Author year]` and resolved in the reference list; each was checked against CrossRef.

## 0. Summary

**What it does.** The package applies one correction algebra — rate correction, baseline subtraction, placebo or vehicle subtraction, a mixed-effects concentration–response model — to hiPSC-CM field-potential data and to clinical 12-lead ECG data, then compares the two on a shared scale: multiples of the free (unbound) clinical Cmax. It will not express an in-vitro endpoint as a QT interval in milliseconds; asking it to raises an error.

**What it shows**, using public data only:

- The clinical side reproduces published analyses: 8 of 8 against [Vicente 2015], 6 of 8 against [Vicente 2019] when stratified by study day, 5 of 6 against [Vicente 2016], and — for the dofetilide arm that had been the project's one unreconciled gap — 2 of 2 to within 0.2 ms against the primary paper [Johannesen 2015]. Our Fridericia QTcF matches the ECGCIPA study's own to a median of 0.000 ms over 5,740 records.
- The in-vitro rate correction reproduces an independent published implementation: 17 of 18 studies within 5 ms of [Millard 2018], median difference 0.65 ms.
- Four drugs are the entire overlap between the in-vitro panel and the single-agent clinical data. At matched free exposure the dish agrees with the clinic in sign for three, overstates the prolongation 1.8- to 3.5-fold, and **inverts** for verapamil (+3.5 ms clinical, −86.3 ms dish).
- The bridge is a property of drugs, not of people: the one subject-matched experiment found no within-person correlation [Blinova 2019].

**What is new** (section 8) is modest and mostly methodological: a symmetric, provenance-typed cascade that refuses the FPD-to-QT conversion; an open, public-data-only side-by-side of dish and clinic; a quantification of informative missingness in the CiPA myocyte data. The J-Tpeak exponent is *not* new. It is published.

**Largest open problems** (section 7): a 5–10 ms between-study shift in the clinical anchor whose remainder is unexplained; informative missingness that published data cannot repair; and a four-drug bridge that cannot support a calibration. The two problems that closed on 4 Nov 2026 — the dofetilide gap and the free-Cmax source — both closed because the user supplied the paywalled papers.

## 1. Question and scope

*Can an in-vitro concentration–response from hiPSC-CM field potentials and a clinical concentration–QTc relationship be passed through the same correction algebra and compared on a common concentration scale — and how far do public data let that be tested?*

Not claimed: that a field potential duration converts into a QT interval; that the dish predicts an individual's response; that either predicts torsade de pointes. The clinical endpoint is QTc prolongation in healthy volunteers, a surrogate, in the population least likely to express the outcome the assay exists to predict.

Context: the Comprehensive in vitro Proarrhythmia Assay (CiPA) paradigm [Sager 2014; Colatsky 2016; Fermini 2016] and its in-silico arm [Li 2019].

## 2. Data

| Source | Content | Size | Notes |
|---|---|---|---|
| CiPA myocyte workbook [Blinova 2018], user-supplied | well-level `ddFPDc`, EAD flag, 28 drugs | 8,859 wells, 15 dataset × platform combinations, 10 sites | Already baseline-, vehicle- and Fridericia-corrected by the source study; the cascade is not re-applied. Validated against all 420 published panels: r = 0.999993 on 1,311 dataset × concentration cells, median absolute difference 0.30 ms. |
| ECGRDVQ [D1], paper [Johannesen 2014b] | dofetilide, quinidine, ranolazine, verapamil, placebo | 22 subjects, crossover | PhysioNet, ODC-By 1.0 |
| ECGDMMLD [D2], paper [Johannesen 2015] | dofetilide, lidocaine+dofetilide, mexiletine+dofetilide, moxifloxacin+diltiazem, placebo | 22 subjects, crossover | PhysioNet, ODC-By 1.0 |
| ECGCIPA [D3], paper [Vicente 2019], design rationale [Vicente 2018] | ranolazine, verapamil, lopinavir/ritonavir, chloroquine (part 1, parallel, 10 per arm); dofetilide ± diltiazem (part 2, crossover, 10 subjects) | 60 subjects | PhysioNet, ODC-By 1.0; ships the study's own derived QTcF, J-Tpeakc and change from baseline |
| Millard 2018 Supplementary Table 1 | baseline beat period, FPD, FPDc for 18 studies | 18 rows | Published as an Enhanced Metafile; recovered by parsing the drawn-text records |
| CiPA Hill fits [W6] | IC50 and Hill coefficient, 12 training drugs × 7 channels | 84 pairs, 19 unfitted | GPL-3.0 with an FDA no-endorsement disclaimer; pinned to a commit |

The clinical anchor is 11 molecules in 104 healthy volunteers, but only four (dofetilide, quinidine, ranolazine, verapamil) are both single-agent and in the in-vitro panel. A prior audit of PhysioNet, Zenodo and the *Scientific Data* descriptors found no other public well-level field-potential dataset.

## 3. Methodology

**3.1 The cascade.** Raw interval → rate-corrected → baseline-subtracted (Δ) → control-subtracted (ΔΔ) → concentration–response model → shared scale. Order is enforced; subtracting a control from an endpoint with no baseline step raises, as does subtracting twice. Records with no matching baseline stay NaN and are counted in the provenance record.

**3.2 Rate correction.** Power family `Ic = I / CL^α` with CL in seconds: Bazett α = 1/2 [Bazett 1920, no DOI], Fridericia α = 1/3 [Fridericia 1920, no DOI], study-specific fitted α, and per-subject α (QTcI) [Malik 2019]. Exponents are fitted on drug-free records only (baseline plus placebo), because fitting on treated records would absorb part of the drug effect. Measured pooled exponents are 0.301 (ECGRDVQ, n = 439, R² = 0.47) and 0.299 (ECGDMMLD, n = 361, R² = 0.40), close to Fridericia's 1/3, which is the expected result in healthy volunteers and a check on the fitting code. Individual exponents were estimable for 22 of 22 and 20 of 22 subjects; the rest fall back to the pooled value, and the fallback is counted.

**3.3 J-Tpeak has two corrections, not interchangeable.** The linear form `JTpc = JTp + 0.150(1 − RR)` is from [Hnatkova 2019]. The power law `JTpc = JTp / (RR/1000)^0.58` is the prespecified formula stated in [Vicente 2019], attributed there to [Johannesen 2014a], which we have not read. We read the exponent from the PMC XML of Vicente 2019, because the plain-text rendering drops superscripts, and confirmed it independently: regressing log(JTp/JTpc) on log(RR) over the 5,740 complete ECGCIPA records recovers 0.5799, intercept 0.000, residual SD 0.15 ms. The same power-law form appears in a comparison of J-Tpeak correction algorithms [Couderc 2017], and the rate dependence of JT interval sections has its own analysis [Hnatkova 2017] (not mined here). The J-Tpeak and Tpeak-Tend intervals come from an automated measurement algorithm described in [Johannesen 2016]. Per record the two corrections differ by a median of 0.44 ms, SD 3.56 ms, maximum 21.1 ms. At the level of drug effects the choice matters little: swapping the power law into ECGRDVQ moves every ΔΔJ-Tpeakc by 1.2 ms or less.

**3.4 Clinical concentration–response.** A linear mixed model [Pinheiro 2000] of ΔΔQTcF on concentration with a subject random intercept, evaluated at the geometric-mean observed Cmax, with the two-sided 90% interval (its upper bound is the one-sided 95% bound that the 10 ms ICH E14 decision uses [ICH E14 Q&A]). Three implementation decisions, each forced by a failure on real data: the predictor is centred at the prediction point, so the effect is the intercept (the delta-method interval had spanned ±5.5×10⁸ ms); the predictor is internally scaled (concentrations span 1.5 to 9,400 ng/mL); and when the between-subject variance is estimated at exactly zero the model is refitted as OLS with cluster-robust errors (4 of 24 clinical fits). The interval uses a *t* quantile on n_groups − 1 degrees of freedom; over 200 simulated crossover studies the coverage of the effect at Cmax was 0.865 (0.845 with a normal quantile) against a nominal 0.90.

**3.5 In-vitro concentration–response.** The same model with the dataset (site × cell line × platform) as the random effect, and **log₁₀(multiples of free Cmax)** as the predictor, because tested concentrations span up to 3,500-fold and a linear fit is dominated by its top point (quinidine returned a negative slope alongside a +129 ms effect at 1×). Missingness is informative (section 6), so every fit is reported twice: over all concentrations and over indices 1–3.

**3.6 The bridge.** ng/mL → nM using the PubChem molecular weight → × fraction unbound → ÷ free Cmax. Both sides are evaluated at the geometric-mean *observed* clinical free Cmax, not a nominal 1×: quinidine reached only 0.29× and verapamil 0.50× of the reference. The geometric mean is taken across subjects of each subject's own maximum post-dose concentration. No in-vitro threshold equivalent to 10 ms exists, so `decision_agreement` requires the caller to supply one and records it as a user convention.

**3.7 Provenance and guardrails.** Every derived value is an `Endpoint` carrying an ordered tuple of `CorrectionRecord`s (step, method, formula, parameters, citation). The modelling layer checks the tuple and names what is missing. `Endpoint.as_qt_milliseconds()` on an in-vitro endpoint raises `DomainError`.

## 4. Code workflow

| Stage | Script | Output |
|---|---|---|
| Build well table and concentration lookup from the workbook | `build_wells.py` | `cipa_myocyte_wells.csv`, `cipa_concentration_lookup.csv` |
| Clinical cascade and concentration–QTc, ECGRDVQ and ECGDMMLD | `clinical_reference.py` | `concentration_qtc_results.csv` |
| ECGCIPA: per-record validation, day-stratified fits, cross-study replication | `ecgcipa_reference.py` | `ecgcipa_concentration_qtc.csv`, `cross_study_replication.csv`, `ecgcipa_correction_validation.csv` |
| In-vitro fits and the dish-versus-clinic table | `invitro_and_bridge.py` | `invitro_concentration_response.csv`, `dish_vs_clinic_bridge.csv` |
| In-vitro rate-correction check against Millard 2018 | `validate_invitro_correction.py` | `invitro_correction_validation.csv` |
| Placebo common-mode diagnostic and J-Tpeak sensitivity | `placebo_commonmode_test.py` | `placebo_commonmode_test.csv`, `ecgrdvq_jtpeakc_power_sensitivity.csv` |

Implementation is in Python with statsmodels [Seabold 2010], SciPy [Virtanen 2020], NumPy [Harris 2020] and pandas [McKinney 2010]; the public clinical data are served by PhysioNet [Goldberger 2000]. The package is layered: `provenance` → `corrections` → `cascade` → `models` → `bridge`, with `io` for PhysioNet, MEA and metafile readers. Tests: 86, collected by `pytest --collect-only` (cascade 10, corrections 14, cross-study diagnostic 5, ECGCIPA io 4, ECGCIPA validation 10, ECGRDVQ validation 14, guardrails 8, in-vitro validation 5, J-Tpeak 6, models 10). Rate corrections are checked against their defining formulas, fitted exponents are recovered from simulated data to 1e-6, the double-delta recovers a known effect through subject offsets and a circadian term, and interval coverage is tested as a frequency over simulated studies.

## 5. Results

**5.1 Validation.**

| Check | Result |
|---|---|
| ECGRDVQ vs [Vicente 2015] Table 2 | 8 of 8 (every published point inside our 90% interval; three "no relationship" cells contain zero) |
| ECGCIPA vs [Vicente 2019], day-stratified | 6 of 8; both misses are verapamil ΔΔQTcF, in opposite directions (5.7 ms low on day 1, 3.3 ms high on day 3), a linear model on a curved response |
| ECGDMMLD vs [Vicente 2016] Table 3 | 5 of 6 |
| ECGDMMLD vs the **primary** paper [Johannesen 2015], dofetilide alone | **2 of 2 to within 0.2 ms**: evening dose 39.1 published vs 39.0 ours; afternoon dose 22.4 vs 22.2. The former "unreconciled 6-7 ms gap" was a disagreement between two published analyses of the same data |
| Per-record QTcF vs the ECGCIPA study's own | median 0.000 ms, SD 0.190 ms; 3 of 5,740 differ by more than 0.05 ms |
| Per-record baseline subtraction vs their change from baseline | median 0.043 ms over 1,915 records |
| In-vitro rate correction vs [Millard 2018] | 17 of 18 within 5 ms, median 0.65 ms; implied exponent 0.3244 vs 0.3333; residual tracks beat-period variability (Spearman ρ = 0.592, p = 0.010), which is Jensen's inequality |

Study 9 (iCell², CardioECR) is 22.0 ms off with a low beat-period coefficient of variation of 0.061, so averaging order does not explain it; the discrepancy is asserted in a test rather than hidden.

**5.2 The bridge (four drugs, the entire overlap).**

| Drug | CiPA risk | × free Cmax | Clinic ΔΔQTcF | Dish ΔΔFPDc | Dish ÷ clinic |
|---|---|---|---|---|---|
| Dofetilide | high | 1.07 | +73.1 [68.5, 77.7] | +132.0 [93.1, 170.8] | 1.81 |
| Quinidine | high | 0.29 | +74.9 [69.3, 80.5] | +135.1 [82.2, 187.9] | 1.80 |
| Ranolazine | low | 0.96 | +11.2 [7.6, 14.8] | +38.8 [22.8, 54.7] | 3.46 |
| Verapamil | low | 0.50 | +3.5 [−0.1, 7.0] | −86.3 [−98.9, −73.6] | −24.7 |

The two pure hERG blockers are inflated by the same factor (1.81, 1.80), which looks like a property of the assay. The dish confidence intervals are 25–106 ms wide against 7–11 ms clinically, because the dominant variance component is between site × platform × cell-line combinations, not between wells.

**5.3 Why verapamil inverts: the channel profile.** Using the CiPA Hill fits at each drug's reference free Cmax, dofetilide and verapamil block hERG by almost the same fraction (23.0% and 22.9%), yet their clinical effects differ twenty-fold. Verapamil also blocks ICaL by 27%, which offsets the hERG block in the clinic; in the dish the calcium arm dominates. Mexiletine (INaL 25%) and ranolazine (INaL 21%) show the same offsetting through late sodium current, and the clinical J-Tpeakc result reproduces it: mexiletine added to dofetilide takes ΔΔQTcF from +37.5 to +20.8 ms and ΔΔJ-Tpeakc from +20.7 to +2.6 ms (each arm at its own dofetilide Cmax). The 28-drug table is in `drug_channel_ecg_profile.csv`; channel data exist for 12 drugs only. For the physics of why a field potential does not resemble an action potential, see [Ernault 2024].

**5.4 Stability of the clinical anchor.** Ranolazine, verapamil and dofetilide were measured in both ECGRDVQ and ECGCIPA. At a matched plasma concentration:

| Drug | Conc. | ECGRDVQ | ECGCIPA | Shift | Shift with no placebo subtraction |
|---|---|---|---|---|---|
| Ranolazine | 2043 ng/mL | 11.2 [7.6, 14.8] | 19.4 [13.1, 25.8] | +8.2 | +6.0 |
| Verapamil | 114 ng/mL | 3.5 [−0.1, 7.0] | 8.4 [6.8, 10.0] | +5.0 | +0.6 |
| Dofetilide | 1.32 ng/mL | 34.1 [29.8, 38.3] | 43.7 [39.0, 48.5] | +9.7 | +7.0 |

ECGCIPA reads higher every time and for dofetilide the intervals do not overlap. ECGCIPA is mixed in design: ranolazine and verapamil are parallel arms, dofetilide a ten-subject crossover with no placebo period of its own, and one part 1 placebo group is subtracted from every arm. The mean post-dose placebo change is −7.6 ms in ECGRDVQ and −9.05 ms in ECGCIPA's single group, so that shared group explains roughly 2–4 ms of the shift, and essentially all of it for verapamil. **For ranolazine and dofetilide about 6–7 ms remains and is unexplained.** An earlier explanation in this project, unpaired placebo subtraction in a parallel design, is withdrawn.

**5.5 The individual-level boundary.** [Blinova 2019] reprogrammed iPSC-CMs from 16 of the volunteers in ECGDMMLD. Subject-specific ΔΔAPD₉₀c slopes did not correlate with the same subjects' clinical ΔΔQTc slopes (Spearman ρ = −0.106, p = 0.70 for dofetilide; ρ = −0.089, p = 0.75 for moxifloxacin), nor did baseline APD₉₀c with baseline QTc (ρ = −0.025). The informative control: the subjects' clinical responses to the two drugs did not correlate with each other (ρ = 0.096, p = 0.73), so the weak signal is on the clinical side as well. Two earlier studies had reported positive subject-level correlations: [Stillitano 2017] in panels built from extreme sotalol responders, and [Shinozawa 2017] in ten healthy volunteers given moxifloxacin. The literature is not settled, and this package therefore carries the result as a boundary: drug-level use only.

## 6. Sources of uncertainty

Magnitudes are from this repository unless marked otherwise; "direction" is the direction of bias where it is known.

| # | Source | Magnitude / evidence | Direction and handling |
|---|---|---|---|
| **Statistical** | | | |
| 1 | Interval coverage | 0.865 vs nominal 0.90 for the effect at Cmax (200 simulated studies); slope 0.900 | Intervals slightly too narrow. Ranolazine's ΔΔJ-Tc call on day 3 of ECGCIPA flips at 10 ms on precision alone (upper bound 9.3 vs 12.0 ms) with point estimates 2 ms apart |
| 2 | Degrees-of-freedom rule | t on n_groups − 1 here; Kenward–Roger [Kenward 1997] is commonly pre-specified in statistical analysis plans [W4] | Not compared |
| 3 | Model simplification | The white-paper model [Garnett 2018] uses ΔQTc with a treatment indicator, nominal-time effects and a centred baseline covariate [W3, W4]; ours fits pre-formed ΔΔ on concentration with a subject intercept | Effect on estimates not quantified. A candidate cause of the unreconciled gaps below |
| 4 | Random-effect collapse | 4 of 24 clinical fits had zero between-subject variance and fall back to OLS with cluster-robust errors | Flagged per fit |
| 5 | Linear model on a curved response | Verapamil ECGCIPA misses: −5.7 ms (day 1), +3.3 ms (day 3) | Emax alternatives exist in the package but are not the default |
| 6 | Not assessed | Residual diagnostics, normality of random effects, influence analysis | Listed in METHODS §7 |
| **Measurement and correction** | | | |
| 7 | Rate-correction choice | Fitted exponents 0.301 and 0.299 vs 1/3; per-record QTcF matches the study's own to 0.19 ms SD | Negligible in healthy volunteers at resting rates; hysteresis not applied (public tables carry one RR per ECG) |
| 8 | J-Tpeak correction form | Per record: median 0.44, SD 3.56, max 21.1 ms. Drug-level: ≤ 1.2 ms | Reported both ways |
| 9 | Triplicate handling | 0.26 ms between mean-of-corrected and corrected-of-mean (Jensen) | Ignored |
| 9b | ECG interval extraction algorithm | **Measured** on PTB-XL+ [Strodthoff 2023], where three algorithms measure the same 21,795 records: the two commercial ones (12SL, Uni-G) differ on QT by a mean of −1.87 ms (SD 7.22) on 9,476 normal records, but 21.3% of all records differ by more than 10 ms. The difference is **rate-dependent**: −0.0196 ms per ms of RR (p = 1e-144), i.e. −5.9 ms across a 300 ms RR change. Fitting the cross-sectional QT/RR exponent on the same records gives 0.3355 (12SL) against 0.3734 (Uni-G); applied to a QT of 400 ms those two exponents put QTc 7.1 ms apart at 90 bpm | A constant bias cancels in ΔΔ; a **rate-dependent** one does not when heart rate moves between baseline and on-drug. Candidate contributor to open problem 1, not a demonstration |
| 10 | Baseline definition | One pre-dose baseline per subject-period, not day-matched across periods | Limitation of the public tables; not quantified |
| **Design and anchor** | | | |
| 11 | Between-study shift | 5.0 to 9.7 ms at matched concentration; 2–4 ms attributable to the shared placebo group, 6–7 ms unexplained | Open problem 1 |
| 12 | ~~Unreconciled dofetilide gap~~ **RESOLVED** | Our 39.0 / 22.2 ms match [Johannesen 2015]'s own 39.1 / 22.4 ms at the evening and afternoon doses. Vicente 2016 Table 3's 45.2 ms disagrees with the primary paper, not with us | No longer an uncertainty in this pipeline; the disagreement between the two published sources is theirs |
| 12b | Random-effect specification | [Johannesen 2015] specifies a random effect on intercept **and** concentration. Refitting with a random slope moves the mean absolute gap over nine published comparisons from 3.16 to 3.04 ms — quinidine improves (-4.0 to -1.0), dofetilide worsens (-7.7 to -11.8) | Tested and rejected as an explanation; intercept-only retained, `random_slope=True` available |
| 13 | Cmax estimate | Geometric mean of each subject's *observed* maximum; sparse sampling biases it low. Ours vs published means: dofetilide 0%, quinidine −2.6%, ranolazine −11%, verapamil −13% | Biased low; test tolerance 15% |
| 14 | Exposure range | Quinidine reached 0.29× and verapamil 0.50× of reference free Cmax | Evaluated at the observed point, not 1× |
| 15 | Population and outcome | Healthy adults aged 19–35, no TdP events in the corpus | Calibrates a surrogate in the population least likely to express the outcome |
| **In-vitro data** | | | |
| 16 | Informative missingness | 2,067 of 8,859 values null; 1,887 are exactly the EAD wells; missingness 0.9% → 4.7% → 26.8% → 60.8% over concentration indices 1–4 | Complete-case slopes biased toward zero. Bounded, not fixed: every fit reported on all concentrations and on indices 1–3 |
| 17 | Consequence for slopes | Quinidine slope −28.3 vs +25.3 ms/decade between the two analyses; terfenadine also flips (2 of 28 drugs). The effect at 0.288× is 135.1 vs 119.4 ms | Quote the effect (stable to about ±16 ms), never the quinidine slope |
| 18 | Source-table defects | 15 under-coded wells at the top concentration for diltiazem, nitrendipine, verapamil (one site, iCell², CardioECR); clozapine / site 1 / iCell² / index 4 has five measurable wells but no plotted value | Documented in the validation report |
| 19 | Mixed endpoints | 2 of 15 datasets (15.2% of wells) are voltage-sensing optical; their `ddFPDc` column holds ΔΔAPD₉₀c | Flagged per well |
| 20 | Site and platform variance | Dish interval widths 25–106 ms vs 7–11 ms clinical | The dominant variance component |
| 20b | Recording-time drift | In one iCell² MEA study, beat rate rose 22.6% and FPD fell 7.7% over a 20-minute equilibration, and wells were excluded for not beating at the desired rate or for low signal [Guerrelli 2024] | Timing of the post-drug read-out differs by site; not modelled |
| 21 | Study 9 | 22.0 ms discrepancy in the rate-correction check, unexplained | Asserted in a test |
| **Translation** | | | |
| 22 | FPD is not QT | Dish overstates 1.8–3.5-fold; verapamil inverts | Structural refusal in code; comparison on concentration, not time |
| 23 | Nominal vs free concentration in the well | **Now mined** [Schocken 2018]: media composition alone moves the baseline. In serum-containing medium iPSC-CMs beat 1.3x slower (interbeat 2102 ± 200 ms vs 1643 ± 143 ms) and baseline FPDc is **575 ± 59 ms against 471 ± 44 ms serum-free** — a 104 ms baseline difference from medium alone. Serum-free medium is used precisely to approximate free drug concentration | **Direction is not what this report previously assumed.** Of 18 drugs that induced arrhythmias, 8 did so at *lower* concentrations in serum-containing medium — the opposite of what protein binding alone predicts, and the paper calls it counterintuitive. Incomplete solubility (quinidine, moxifloxacin) also makes nominal exceed actual. The earlier claim that correcting for free concentration would simply widen the dish-clinic discrepancy is **withdrawn**; the sign is not established |
| 24 | Unbound fraction | Single literature point estimates; four CiPA drugs have none because [Crumb 2016] and [Li 2017] could not be retrieved; fu changes with acute-phase protein levels | Blank, never guessed |
| 25 | Reference Cmax | [Blinova 2018] and the CiPA training table [Dutta 2017] disagree: terfenadine 13.3×, mexiletine 1.65×, verapamil 1.62×, ondansetron 0.38× | Bridge unaffected (both sides use one reference); channel-block percentages are affected |
| 26 | Immaturity and scale | Low IK1, spontaneous beating, no autonomic input, no conduction system | Why FPD and QT are analogous and not convertible |
| 27 | Individual prediction | No within-person correlation [Blinova 2019]; positive earlier reports [Stillitano 2017; Shinozawa 2017] | Drug-level use only |

## 7. Open problems

1. **The 6–7 ms between-study shift** for ranolazine and dofetilide. *Would settle it:* refit both studies with the full white-paper model (uncertainty 3); refit ECGCIPA day by day; compare ECG extraction methods and population composition; ideally obtain waveforms. *Partial progress:* uncertainty 9b now shows that interval-extraction method alone carries a rate-dependent QT bias of about 6 ms across a 300 ms RR change, and changes the fitted QT/RR exponent from 0.336 to 0.373 on identical records. That is the right order of magnitude, so extraction cannot be dismissed — but it is not evidence that it *is* the cause, because the extraction method used by each clinical study was not established here, and the drugs changed mean heart rate by less than 10 bpm [Vicente 2019], which spans only about a 140 ms RR change.
2. ~~**The ECGDMMLD dofetilide gap**~~ **CLOSED.** [Johannesen 2015] was obtained and reports 39.1 ms (evening) and 22.4 ms (afternoon) for dofetilide alone; this pipeline gives 39.0 and 22.2. The gap was against Vicente 2016 Table 3, which disagrees with its own study's primary paper. The random-slope specification that paper prescribes was tested and does not explain the difference (uncertainty 12b).
3. **Informative missingness cannot be repaired.** The underlying durations do not exist in the published data. *Next step:* a tipping-point analysis that assigns the missing wells an FPD of the observed maximum plus δ and reports the δ at which a conclusion changes. Not yet done.
4. **A four-drug bridge cannot support a calibration.** *Would settle it:* more single-agent clinical anchors with open concentration–QTc data. None are public; the 34-drug retrospective [Johannesen 2014a] is a summary analysis.
5. **Dish-side free concentration.** Partly addressed: [Schocken 2018] shows medium composition moves baseline FPDc by ~104 ms and that serum's effect on drug sensitivity runs counter to a pure protein-binding argument for 8 of 18 drugs. What remains is that free concentration in the well was still not measured per drug per site, so the dish axis remains nominal. *Would settle it:* per-well free-concentration measurement.
6. ~~**Unbound fractions and reference Cmax**~~ **Largely closed.** [Crumb 2016] Table 1 supplies free Cmax in nM for all 30 drugs of its panel (`data/crumb2016_free_cmax.csv`), which is the quantity the bridge actually consumes; the Blinova reference matches it for 11 of 12 shared drugs (median ratio 1.000, quinidine the exception). The bridge is in any case invariant to the choice, since the same reference scales both axes. What is still missing is the *fraction unbound* itself for bepridil, chlorpromazine, cisapride and terfenadine — neither Crumb nor [Li 2017] tabulates fu, only free Cmax.
7. **Model form:** curvature (verapamil), an effect compartment for hysteresis, and RR hysteresis. *Would settle it:* day-stratified Emax fits and beat-level RR.
8. **Interval validity:** the degrees-of-freedom rule, the 0.865 coverage shortfall, residual diagnostics.
9. **Validating the in-vitro baseline and vehicle steps.** No public uncorrected well-level dataset exists, so these are tested only on synthetic data and, for the shared code path, on the clinical corpus.
10. **The individual-level literature conflict** between [Stillitano 2017], [Shinozawa 2017] and [Blinova 2019] is unresolved. This package does not need to resolve it, only to stay on the drug-level side of it.
11. **From effect to decision.** No in-vitro threshold equivalent to the clinical 10 ms exists, and no arrhythmic outcome exists in the corpus.

## 8. What is novel for us

**Method of the screen.** Fourteen short PubMed queries, seven OpenAlex queries and three web queries (November 2026), plus a full-text check of [Vicente 2019] for the J-Tpeak exponent. A "not found" here is a statement about those sources only. Earlier screens in this programme were overturned by sources the first pass missed, so none of the verdicts below is a claim of priority.

| Claim | Verdict | Closest prior work | Note |
|---|---|---|---|
| A symmetric correction cascade across hiPSC-CM and clinical data, with typed provenance and a structural refusal to express FPD as QT | **Not found elsewhere; a software-design contribution** | Open MEA tools cover signal processing only [Dunham 2022; Lee 2024]; the open concentration–QTc code found is an R package driving the commercial MonolixSuite [W1] and an R power-calculation repository [W2] | The README's earlier claim that concentration–QTc analysis has "no open implementation" is withdrawn |
| An open, public-data-only side-by-side of dish and clinic at free-Cmax multiples, with the clinical side re-derived from open ECG data | **Not found; modest** | Categorical comparison [Blinova 2016]; subject-level comparison [Blinova 2019]; a high-throughput hiPSC-CM assay built to predict clinical ECG changes [Kilfoil 2021] (title only; not read); a MEA study that estimated repolarisation at clinical Cmax from the dose–response curve [W5] | Effect-at-Cmax in the dish is not new. The new part is the symmetric re-derivation of the clinical side |
| Quantified informative missingness in the CiPA myocyte data (1,887 of 2,067 nulls are EAD wells; 0.9% → 60.8%) and the double-fit bound | **Not found; a quantification, not a discovery** | That arrhythmic wells lack a measurable FPD is implicit in the design; how [Blinova 2018] handled it was not checked line by line | Theory of the problem: [Little 2019] |
| The J-Tpeakc exponent of 0.58 | **Not novel; now confirmed at source** | [Johannesen 2014a], read 4 Nov 2026: "the coefficient derived from data in this study: J-Tpeakc = J-Tpeak / RR^0.58", pooled over 34 thorough-QT studies | Ours: a verification on 5,740 records, the 21 ms per-record spread between the two forms, and the ≤ 1.2 ms drug-level effect |
| The ECGRDVQ–ECGCIPA shift at matched concentration, with a placebo decomposition | **Not found; an empirical observation with an unexplained remainder** | Design and baseline effects on concentration–QTc [Garnett 2016] | The remainder is open problem 1 |
| The individual-level null | **Not ours** | [Blinova 2019]; positive earlier reports [Stillitano 2017; Shinozawa 2017] | Ours: carrying it as a hard boundary on the package |
| Recovering a published table from an Enhanced Metafile; numerical traps in the mixed model | **Useful engineering, not scientific novelty** | | |

## 8.1 Appendix: other public ECG corpora, and why they do not widen the bridge

Eight widely used ECG datasets were assessed. The binding requirement for this project is **paired drug exposure and ECG** — a measured plasma concentration beside an interval, with a baseline and a control arm. None of these has any of that: a keyword scan of their PhysioNet pages for *drug*, *dose*, *placebo*, *concentration*, *plasma* and *pharmacokinetic* returned nothing (PTB has a single mention of "medication"). They are diagnosis and arrhythmia-classification corpora, built for a different task.

| Corpus | What it is | Why it cannot widen the bridge | Residual use here |
|---|---|---|---|
| MIT-BIH Arrhythmia [Moody 2001] | 48 half-hour ambulatory records, 47 subjects, 2 channels, 360 Hz, ~110,000 beat annotations | 2-channel Holter, no 12-lead QT measurement, no drug data | Long records give within-subject rate range; could support hysteresis work, which the clinical tables cannot |
| PhysioNet/CinC 2017 [Clifford 2017] | short single-lead recordings, AF classification | single lead; no QT, no drug | None |
| INCART [D4] | 75 records of 30 min, 12-lead, 257 Hz, from 32 Holter recordings; coronary-artery-disease patients | no drug exposure; 257 Hz is coarse for interval onsets | Same hysteresis use as MIT-BIH |
| CPSC-2018 [Liu 2018] | 12-lead clinical ECGs, 9 diagnostic classes | diagnosis labels only, no drug | None |
| PTB [Bousseljot 1995] | 549 records, 290 subjects, 15 leads at 1000 Hz | no drug exposure; small | High sampling rate; waveform work only |
| PTB-XL [Wagner 2020] | 21,799 clinical 12-lead ECGs, 18,869 patients, 10 s, SCP statements | 10 s per patient, cross-sectional, no drug | Normal-population reference; too short for within-subject QT/RR |
| PhysioNet/CinC 2020-2021 [Perez Alday 2020; Reyna 2021] | large multi-source 12-lead collection, multi-label diagnosis | aggregates the above; no drug | None |
| CCDD [Zhang 2010] | large Chinese 12-lead diagnostic database | diagnosis labels; no drug; access is not open on PhysioNet | None |

**The one exception is PTB-XL+** [Strodthoff 2023], which supplements PTB-XL with ECG features from two commercial algorithms (12SL, Uni-G) and one open-source one (ECGDeli [Pilia 2021]), harmonised to common names and units, **for the same records**. That makes it a controlled experiment on measurement method rather than on drugs, and it answers a question this report had left open. The result is uncertainty 9b: a small constant QT bias between the two commercial algorithms but a rate-dependent one of about 6 ms across a 300 ms RR change, and a fitted QT/RR exponent of 0.3355 against 0.3734 on identical records. ECGDeli's mapped QT feature is a maximum across leads (`QTi_max`) and reads a median 108 ms longer, so it is not comparable and was excluded rather than treated as a third estimate.

Two limits on that measurement. PTB-XL is patients rather than healthy volunteers, and each record is a single 10-second resting ECG, so the exponent fitted here is **cross-sectional** — it mixes between-subject differences with rate dependence and is not the within-subject exponent that `individual_correction` estimates. And the comparison speaks only to these two algorithms, not to the specific software used by the clinical studies, which was not established.

The general conclusion: the public ECG landscape is large but almost entirely built for diagnosis from waveforms. Concentration-paired ECG remains the scarce resource, and the four-drug overlap in open problem 4 is not relieved by any of these.

## 9. Limits of this report

- The repository is not on GitHub; no credential is configured.
- All seven previously closed or unmined papers were supplied by the user on 4 Nov 2026 and have now been read: [Johannesen 2014a], [Johannesen 2014b], [Johannesen 2015], [Crumb 2016], [Li 2017], [Strauss 2019], [Schocken 2018]. [Li 2017] and [Strauss 2019] were read for the specific values sought (fraction unbound; CiPA status) and not mined in full.
- The PTB-XL+ measurement is observational and cross-sectional; it bounds what extraction method *can* do, not what it *did* do in the clinical studies.
- Channel data exist for 12 of the 28 panel drugs. For the other 16 the characterisation in `drug_channel_ecg_profile.csv` is general pharmacology and was not retrieved.
- Corrected during this project, so the record is visible: a J-Tpeak formula credited to the wrong first author; a cross-study explanation that was only partly right and a design label that was wrong; a "discovered" exponent that was published; a "no open implementation" claim; a test count taken from a hand count; and a quinidine caveat that named the wrong weakness.

## References

- **[Blinova 2016]** Blinova K et al. (2016) Comprehensive Translational Assessment of Human-Induced Pluripotent Stem Cell Derived Cardiomyocytes for Evaluating Drug-Induced Arrhythmias. *Toxicological Sciences* 155:234-247. doi:10.1093/toxsci/kfw200
- **[Blinova 2018]** Blinova K et al. (2018) International Multisite Study of Human-Induced Pluripotent Stem Cell-Derived Cardiomyocytes for Drug Proarrhythmic Potential Assessment. *Cell Reports* 24:3582-3592. doi:10.1016/j.celrep.2018.08.079
- **[Blinova 2019]** Blinova K et al. (2019) Clinical Trial in a Dish: Personalized Stem Cell-Derived Cardiomyocyte Assay Compared With Clinical Trial Results for Two QT-Prolonging Drugs. *Clinical and Translational Science* 12:687-697. doi:10.1111/cts.12674
- **[Bousseljot 1995]** Bousseljot R et al. (1995) Nutzung der EKG-Signaldatenbank CARDIODAT der PTB über das Internet. *Biomedizinische Technik/Biomedical Engineering*:317-318. doi:10.1515/bmte.1995.40.s1.317
- **[Clifford 2017]** Clifford G et al. (2017) AF Classification from a Short Single Lead ECG Recording: the Physionet Computing in Cardiology Challenge 2017. *Computing in Cardiology Conference (CinC)* 44. doi:10.22489/cinc.2017.065-469
- **[Colatsky 2016]** Colatsky T et al. (2016) The Comprehensive in Vitro Proarrhythmia Assay (CiPA) initiative — Update on progress. *Journal of Pharmacological and Toxicological Methods* 81:15-20. doi:10.1016/j.vascn.2016.06.002
- **[Couderc 2017]** Couderc JP et al. (2017) An evaluation of multiple algorithms for the measurement of the heart rate corrected JTpeak interval. *Journal of Electrocardiology* 50:769-775. doi:10.1016/j.jelectrocard.2017.08.025
- **[Crumb 2016]** Crumb WJ et al. (2016) An evaluation of 30 clinical drugs against the comprehensive in vitro proarrhythmia assay (CiPA) proposed ion channel panel. *Journal of Pharmacological and Toxicological Methods* 81:251-262. doi:10.1016/j.vascn.2016.03.009
- **[Dunham 2022]** Dunham CS et al. (2022) Cardio PyMEA: A user-friendly, open-source Python application for cardiomyocyte microelectrode array analysis. *PLOS ONE* 17:e0266647. doi:10.1371/journal.pone.0266647
- **[Dutta 2017]** Dutta S et al. (2017) Optimization of an In silico Cardiac Cell Model for Proarrhythmia Risk Assessment. *Frontiers in Physiology* 8. doi:10.3389/fphys.2017.00616
- **[Ernault 2024]** Ernault AC et al. (2024) Interpretation of field and LEAP potentials recorded from cardiomyocyte monolayers. *American Journal of Physiology-Heart and Circulatory Physiology* 326:H800-H811. doi:10.1152/ajpheart.00463.2023
- **[Fermini 2016]** Fermini B et al. (2016) A New Perspective in the Field of Cardiac Safety Testing through the Comprehensive In Vitro Proarrhythmia Assay Paradigm. *SLAS Discovery* 21:1-11. doi:10.1177/1087057115594589
- **[Garnett 2016]** Garnett C et al. (2016) Operational Characteristics of Linear Concentration-QT Models for Assessing QTc Interval in the Thorough QT and Phase I Clinical Studies. *Clinical Pharmacology & Therapeutics* 100:170-178. doi:10.1002/cpt.361
- **[Garnett 2018]** Garnett C et al. (2018) Scientific white paper on concentration-QTc modeling. *Journal of Pharmacokinetics and Pharmacodynamics* 45:383-397. doi:10.1007/s10928-017-9558-5
- **[Goldberger 2000]** Goldberger AL et al. (2000) PhysioBank, PhysioToolkit, and PhysioNet. *Circulation* 101. doi:10.1161/01.cir.101.23.e215
- **[Guerrelli 2024]** Guerrelli D et al. (2024) hiPSC-CM electrophysiology: impact of temporal changes and study parameters on experimental reproducibility. *American Journal of Physiology-Heart and Circulatory Physiology* 327:H12-H27. doi:10.1152/ajpheart.00631.2023
- **[Harris 2020]** Harris CR et al. (2020) Array programming with NumPy. *Nature* 585:357-362. doi:10.1038/s41586-020-2649-2
- **[Hnatkova 2017]** Hnatkova K et al. (2017) Heart rate dependency of JT interval sections. *Journal of Electrocardiology* 50:814-824. doi:10.1016/j.jelectrocard.2017.08.005
- **[Hnatkova 2019]** Hnatkova K et al. (2019) Heart Rate Correction of the J-to-Tpeak Interval. *Scientific Reports* 9. doi:10.1038/s41598-019-51491-4
- **[Johannesen 2014a]** Johannesen L et al. (2014) Improving the Assessment of Heart Toxicity for All New Drugs Through Translational Regulatory Science. *Clinical Pharmacology & Therapeutics* 95:501-508. doi:10.1038/clpt.2013.238
- **[Johannesen 2014b]** Johannesen L et al. (2014) Differentiating Drug-Induced Multichannel Block on the Electrocardiogram: Randomized Study of Dofetilide, Quinidine, Ranolazine, and Verapamil. *Clinical Pharmacology & Therapeutics* 96:549-558. doi:10.1038/clpt.2014.155
- **[Johannesen 2015]** Johannesen L et al. (2015) Late sodium current block for drug-induced long QT syndrome: Results from a prospective clinical trial. *Clinical Pharmacology & Therapeutics* 99:214-223. doi:10.1002/cpt.205
- **[Johannesen 2016]** Johannesen L et al. (2016) Automated Algorithm for J-Tpeak and Tpeak-Tend Assessment of Drug-Induced Proarrhythmia Risk. *PLOS ONE* 11:e0166925. doi:10.1371/journal.pone.0166925
- **[Kenward 1997]** Kenward MG et al. (1997) Small Sample Inference for Fixed Effects from Restricted Maximum Likelihood. *Biometrics* 53:983. doi:10.2307/2533558
- **[Kilfoil 2021]** Kilfoil P et al. (2021) Characterization of a high throughput human stem cell cardiomyocyte assay to predict drug-induced changes in clinical electrocardiogram parameters. *European Journal of Pharmacology* 912:174584. doi:10.1016/j.ejphar.2021.174584
- **[Lee 2024]** Lee J et al. (2024) CardioMEA: comprehensive data analysis platform for studying cardiac diseases and drug responses. *Frontiers in Physiology* 15. doi:10.3389/fphys.2024.1472126
- **[Li 2017]** Li Z et al. (2017) Improving the In Silico Assessment of Proarrhythmia Risk by Combining hERG (Human Ether-à-go-go-Related Gene) Channel-Drug Binding Kinetics and Multichannel Pharmacology. *Circulation: Arrhythmia and Electrophysiology* 10. doi:10.1161/CIRCEP.116.004628
- **[Li 2019]** Li Z et al. (2019) Assessment of an In Silico Mechanistic Model for Proarrhythmia Risk Prediction Under the CiPA Initiative. *Clinical Pharmacology & Therapeutics* 105:466-475. doi:10.1002/cpt.1184
- **[Little 2019]** Little R et al. (2019) Statistical Analysis with Missing Data, Third Edition. *Wiley Series in Probability and Statistics*. doi:10.1002/9781119482260
- **[Liu 2018]** Liu F et al. (2018) An Open Access Database for Evaluating the Algorithms of Electrocardiogram Rhythm and Morphology Abnormality Detection. *Journal of Medical Imaging and Health Informatics* 8:1368-1373. doi:10.1166/jmihi.2018.2442
- **[Malik 2019]** Malik M et al. (2019) Implications of Individual QT/RR Profiles—Part 1: Inaccuracies and Problems of Population-Specific QT/Heart Rate Corrections. *Drug Safety* 42:401-414. doi:10.1007/s40264-018-0736-1
- **[McKinney 2010]** McKinney W (2010) Data Structures for Statistical Computing in Python. *Proceedings of the Python in Science Conference*:56-61. doi:10.25080/majora-92bf1922-00a
- **[Millard 2018]** Millard D et al. (2018) Cross-Site Reliability of Human Induced Pluripotent stem cell-derived Cardiomyocyte Based Safety Assays Using Microelectrode Arrays: Results from a Blinded CiPA Pilot Study. *Toxicological Sciences* 164:550-562. doi:10.1093/toxsci/kfy110
- **[Moody 2001]** Moody GB et al. (2001) The impact of the MIT-BIH Arrhythmia Database. *IEEE Engineering in Medicine and Biology Magazine* 20:45-50. doi:10.1109/51.932724
- **[Perez Alday 2020]** Perez Alday EA et al. (2020) Classification of 12-lead ECGs: the PhysioNet/Computing in Cardiology Challenge 2020. *Physiological Measurement* 41:124003. doi:10.1088/1361-6579/abc960
- **[Pilia 2021]** Pilia N et al. (2021) ECGdeli - An open source ECG delineation toolbox for MATLAB. *SoftwareX* 13:100639. doi:10.1016/j.softx.2020.100639
- **[Pinheiro 2000]** Pinheiro JC et al. (2000) Mixed-Effects Models in S and S-PLUS. *Statistics and Computing*. doi:10.1007/978-1-4419-0318-1
- **[Reyna 2021]** Reyna MA et al. (2021) Will Two Do? Varying Dimensions in Electrocardiography: The PhysioNet/Computing in Cardiology Challenge 2021. *2021 Computing in Cardiology (CinC)*:1-4. doi:10.23919/cinc53138.2021.9662687
- **[Sager 2014]** Sager PT et al. (2014) Rechanneling the cardiac proarrhythmia safety paradigm: A meeting report from the Cardiac Safety Research Consortium. *American Heart Journal* 167:292-300. doi:10.1016/j.ahj.2013.11.004
- **[Schocken 2018]** Schocken D et al. (2018) Comparative analysis of media effects on human induced pluripotent stem cell-derived cardiomyocytes in proarrhythmia risk assessment. *Journal of Pharmacological and Toxicological Methods* 90:39-47. doi:10.1016/j.vascn.2017.11.002
- **[Seabold 2010]** Seabold S et al. (2010) Statsmodels: Econometric and Statistical Modeling with Python. *Proceedings of the Python in Science Conference*:92-96. doi:10.25080/majora-92bf1922-011
- **[Shinozawa 2017]** Shinozawa T et al. (2017) Recapitulation of Clinical Individual Susceptibility to Drug-Induced QT Prolongation in Healthy Subjects Using iPSC-Derived Cardiomyocytes. *Stem Cell Reports* 8:226-234. doi:10.1016/j.stemcr.2016.12.014
- **[Stillitano 2017]** Stillitano F et al. (2017) Modeling susceptibility to drug-induced long QT with a panel of subject-specific induced pluripotent stem cells. *eLife* 6. doi:10.7554/elife.19406
- **[Strauss 2019]** Strauss DG et al. (2019) Comprehensive In Vitro Proarrhythmia Assay (CiPA) Update from a Cardiac Safety Research Consortium / Health and Environmental Sciences Institute / FDA Meeting. *Therapeutic Innovation & Regulatory Science* 53:519-525. doi:10.1177/2168479018795117
- **[Strodthoff 2023]** Strodthoff N et al. (2023) PTB-XL+, a comprehensive electrocardiographic feature dataset. *Scientific Data* 10. doi:10.1038/s41597-023-02153-8
- **[Vicente 2015]** Vicente J et al. (2015) Comprehensive T wave Morphology Assessment in a Randomized Clinical Study of Dofetilide, Quinidine, Ranolazine, and Verapamil. *Journal of the American Heart Association* 4. doi:10.1161/JAHA.114.001615
- **[Vicente 2016]** Vicente J et al. (2016) Electrocardiographic Biomarkers for Detection of Drug-Induced Late Sodium Current Block. *PLOS ONE* 11:e0163619. doi:10.1371/journal.pone.0163619
- **[Vicente 2018]** Vicente J et al. (2018) Mechanistic Model-Informed Proarrhythmic Risk Assessment of Drugs: Review of the “CiPA” Initiative and Design of a Prospective Clinical Validation Study. *Clinical Pharmacology & Therapeutics* 103:54-66. doi:10.1002/cpt.896
- **[Vicente 2019]** Vicente J et al. (2019) Assessment of Multi-Ion Channel Block in a Phase I Randomized Study Design: Results of the CiPA Phase I ECG Biomarker Validation Study. *Clinical Pharmacology & Therapeutics* 105:943-953. doi:10.1002/cpt.1303
- **[Virtanen 2020]** Virtanen P et al. (2020) SciPy 1.0: fundamental algorithms for scientific computing in Python. *Nature Methods* 17:261-272. doi:10.1038/s41592-019-0686-2
- **[Wagner 2020]** Wagner P et al. (2020) PTB-XL, a large publicly available electrocardiography dataset. *Scientific Data* 7. doi:10.1038/s41597-020-0495-6
- **[Zhang 2010]** Zhang Jw et al. (2010) Chinese Cardiovascular Disease Database (CCDD) and Its Management Tool. *2010 IEEE International Conference on BioInformatics and BioEngineering*:66-72. doi:10.1109/bibe.2010.19

**Non-DOI sources** (identifier absent where stated)

- **[Bazett 1920]** Bazett HC. *Heart* 1920;7:353-370. No DOI.
- **[Fridericia 1920]** Fridericia LS. *Acta Med Scand* 1920;53:469-486. No DOI.
- **[ICH E14 Q&A]** International Council for Harmonisation. E14 Questions & Answers (R3). No DOI; year not verified here.
- **[D1]** ECG Effects of Ranolazine, Dofetilide, Verapamil, and Quinidine v1.0.0. PhysioNet. doi:10.13026/C2HP45. ODC-By 1.0.
- **[D2]** ECG Effects of Dofetilide, Moxifloxacin, Dofetilide+Mexiletine, Dofetilide+Lidocaine and Moxifloxacin+Diltiazem v1.0.0 (title as truncated by the page). PhysioNet. doi:10.13026/C2D016. ODC-By 1.0.
- **[D3]** CiPA ECG Validation Study v1.0.0. PhysioNet. doi:10.13026/C2967M. ODC-By 1.0.
- **[D4]** St Petersburg INCART 12-lead Arrhythmia Database v1.0.0. PhysioNet. doi:10.13026/C2V88N. No associated paper found.
- **[W1]** Simulations Plus / PAGE meeting abstract: Automated concentration-QT data preparation, model selection and reporting in R. https://www.page-meeting.org/Abstracts/automated-concentration-qt-data-preparation-model-selection-and-reporting-in-r/ (authors and year not retrieved).
- **[W2]** cqtcpower: power and sample size for concentration-QTc analyses of parallel-group single ascending dose studies. https://github.com/YINGTIAN-HU/cqtcpower
- **[W3]** Beyond the linear model in concentration-QT analysis. PMC12065759 (journal, authors and year not retrieved).
- **[W4]** Statistical analysis plan, ClinicalTrials.gov NCT04821063, section 10.3.1 (an example of Kenward–Roger degrees of freedom and a treatment-and-baseline model). https://cdn.clinicaltrials.gov/large-docs/63/NCT04821063/SAP_001.pdf
- **[W5]** Assessment of Drug Proarrhythmic Potential in Electrically Paced Human Induced Pluripotent Stem Cell-Derived Ventricular Cardiomyocytes Using Multielectrode Array. https://www.sciencedirect.com/science/article/pii/S2472555222066825 (authors and year not retrieved).
- **[W6]** FDA CiPA repository, ion-channel Hill fits, commit 539ff80d6616895de5d943e7bfb77b346d195db4. GPL-3.0 with an FDA no-endorsement disclaimer.
