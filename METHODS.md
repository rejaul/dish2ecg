# Methods

Every equation the package applies, the source it comes from, and the decision
behind it. Notation: `I` is the interval being corrected, `CL` the cycle length
governing it. In the clinic `I` = QT or J-Tpeak and `CL` = RR; in the dish `I` =
FPD or APD₉₀ and `CL` = the interbeat interval. Inputs and outputs are in
milliseconds; the classical formulas are defined with the cycle length in
seconds, and the conversion is internal and recorded in the provenance string.

## 1. Rate correction

### 1.1 Power family

    Ic = I / (CL_s)^α

| α | name | source |
|---|---|---|
| 1/2 | Bazett | Bazett HC, *Heart* 1920;7:353-370 |
| 1/3 | Fridericia | Fridericia LS, *Acta Med Scand* 1920;53:469-486; the ICH E14 default |
| fitted, pooled | study-specific | log-linear regression of log I on log CL over drug-free records |
| fitted, per subject | individual (QTcI) | Malik M et al., *Drug Saf* 2019;42:401-414 |

Bazett is retained for comparability with older reports and over-corrects above
roughly 70 bpm; the package computes it but does not recommend it.

Fitted exponents are estimated on **drug-free records only** — baseline records
plus the whole placebo arm — which is how study-specific and population
corrections are defined. Fitting them on on-drug records would absorb part of the
drug effect into the correction.

Measured on the two public studies, the pooled exponent is 0.301 (ECGRDVQ,
n = 439 drug-free records, R² = 0.47) and 0.299 (ECGDMMLD, n = 361, R² = 0.40) —
both close to Fridericia's 1/3, which is the expected result in healthy
volunteers and a check that the fitting code is doing what it claims. Individual
exponents were estimable for 22 of 22 subjects in ECGRDVQ and 20 of 22 in
ECGDMMLD (minimum 10 usable records); subjects below the threshold fall back to
the pooled exponent and the fallback is counted in the provenance record, never
applied silently.

### 1.2 Linear family

    Framingham:  QTc = QT + 154 · (1 − RR_s)      Sagie A et al., Am J Cardiol 1992;70:797-801
    Hodges:      QTc = QT + 1.75 · (HR − 60)      Hodges M et al., J Electrocardiol 1983;16:17-24

### 1.3 J-Tpeak

    JTpc = JTp + 0.150 · (1 − RR)          all quantities in seconds

Hnatkova K, Vicente J, Johannesen L, Garnett C, Strauss DG, Stockbridge N,
Malik M, *Sci Rep* 2019;9:15060, doi:10.1038/s41598-019-51491-4,
derived from QT/RR, JTp/RR and JT50/RR profiles in 523 healthy subjects. The
constant was read from the source paper, not recalled. Applied here in
milliseconds as `+150 · (1 − RR_ms/1000)`.

**There are two J-Tpeak corrections in this package and they are not
interchangeable.** The linear form above is Hnatkova 2019 (Sci Rep 9:15060, doi:10.1038/s41598-019-51491-4; Vicente is the second author, so citing it as "Vicente 2019" collides with the unrelated Vicente 2019 CPT ECGCIPA paper). The FDA's CiPA ECG work
uses a power law instead,

    JTpc = JTp / (RR_s)^0.58

The exponent is published, and the source has now been read. Johannesen L et al.,
*Clin Pharmacol Ther* 2014;95:501-508, doi:10.1038/clpt.2013.238, states it directly:
"heart-rate correction was performed only for J-Tpeak (using the coefficient derived
from data in this study: J-Tpeakc = J-Tpeak / RR^0.58) and QT (using Fridericia's
correction)". The coefficient was derived there by pooling subjects from 34 thorough-QT
studies, and the companion trial paper (doi:10.1038/clpt.2014.155) repeats it as
"J-Tpeakc = J-Tpeak/RR^0.58 with RR in seconds". Vicente 2019 restates the same formula
with RR in ms. It was read from the PMC XML
of Vicente 2019, because the plain-text rendering drops superscripts. As an
independent check, PhysioNet
ECGCIPA ships the study's own derived J-Tpeakc next to the raw J-Tpeak and RR;
regressing log(JTp/JTpc) on log(RR) over its 5740 complete records gives a pure
power law with exponent **0.5799**, intercept 0.000 and residual SD 0.15 ms.
`jtpeakc_power` uses 0.58 by default and records that provenance. Over the same
records the two corrections differ by a median of 0.44 ms, SD 3.56 ms and a
maximum of 21.1 ms, which is large against a 10 ms decision threshold: choose
deliberately, and say which you used.

The ECGRDVQ and ECGDMMLD reference script uses the linear form, to follow the
published tables; the ECGCIPA script uses the power law. Swapping the power law
into ECGRDVQ moves every ΔΔJ-Tpeakc by 1.2 ms or less (dofetilide 38.7 → 39.4,
quinidine 27.7 → 27.6, ranolazine 2.8 → 2.4, verapamil 0.6 → −0.6; verapamil's
interval contains zero under both), so the choice changes no conclusion drawn from
those studies even though individual records differ by up to 21 ms.
(`results/ecgrdvq_jtpeakc_power_sensitivity.csv`.)

Two caveats carried in the code. First, the published derivation uses
**hysteresis-corrected** RR; the PhysioNet tables publish one summary RR per
extracted ECG, so the correction as applied here uses the instantaneous RR and is
an approximation. `hysteresis_corrected_rr` is implemented for data that do carry
beat-to-beat history. Second, `QT = QRS + JTPEAK + TPEAKTEND` holds for 100 % of
complete rows in both studies — QT runs from QRS onset and J-Tpeak from QRS
offset — so J-Tpeak is already the QRS-subtracted interval and is reconstructable
where missing.

## 2. The subtraction cascade

    Δ   = corrected value  −  baseline of the same subject (or well) in the same arm
    ΔΔ  = Δ  −  Δ of the control arm        (placebo in the clinic, vehicle in the dish)

Both steps run through one code path, differing only in column names and the
identity of the control arm. Order is enforced: subtracting a control arm from an
endpoint that has not had a baseline removed raises, and so does subtracting
either twice.

**Baseline definition for the public studies.** `BASELINE == 'Y'` is an exact
alias for `TPT == -0.5`, and there is one pre-dose baseline per subject per
period — not a day-matched baseline across periods. The cascade is therefore run
with `match_on_time=False`, taking the subject-by-period baseline mean. Records
with no matching baseline are left as NaN and **counted** in the provenance
record; they are never dropped silently.

**Replicates.** Three ECGs are extracted at every nominal timepoint. They are
averaged per subject × visit × timepoint before modelling, per the ICH E14
triplicate convention (3147 of 3148 cells have 3 replicates; one ECGDMMLD cell
has 2).

## 3. Concentration–response

### 3.1 Clinical

Linear mixed-effects model per ICH E14 Q&A (R3):

    ΔΔQTcF_ij = (β₀ + b_i) + β₁ · C_ij + ε_ij,     b_i ~ N(0, σ²_b)

with subject as the random intercept. Reported: the slope, the predicted effect
at the geometric-mean observed Cmax, and the two-sided 90 % confidence interval
of that prediction. The regulatory question is whether the **upper** bound
exceeds 10 ms.

Three implementation decisions, each made because the naive version failed on
real data:

1. **Centred at the prediction point.** The model is fitted on
   `C − C_target`, so the predicted effect *is* the intercept and its standard
   error is read directly rather than assembled by the delta method from three
   covariance entries. With ~20 subjects the off-diagonal covariance terms are
   unreliable, and the delta-method interval was observed to span ±5.5 × 10⁸ ms
   while the point estimate stayed correct.
2. **Internally scaled predictor.** Concentrations across the panel span 1.5
   ng/mL (dofetilide) to 9400 ng/mL (moxifloxacin). The predictor is unit-scaled
   for fitting and the slope mapped back.
3. **Boundary fallback.** When the between-subject variance is estimated at
   exactly zero — routine after an individual rate correction, which removes most
   between-subject variation — the mixed model has collapsed to OLS, and
   statsmodels returns a degenerate intercept. The package detects this, refits as
   OLS with cluster-robust standard errors by subject, records
   `boundary_fallback=True` and says so in `random_effects`. Four of the 24
   clinical fits take this path.

**Interval quantile.** Two-sided 90 % from a *t* distribution on `n_groups − 1`
degrees of freedom. This is the regulatory convention and it matters at this
sample size: measured over 200 simulated crossover studies, coverage of the true
effect at Cmax was 0.845 with a normal quantile against a nominal 0.90, and 0.865
with the t quantile; slope coverage was 0.900 either way. The residual shortfall
is documented rather than tuned away.

### 3.2 In vitro

Same model, with two changes.

**Random effects** become the dataset — site × cell line × platform — rather than
the subject. Fifteen such datasets exist in the CiPA myocyte study.

**The predictor is log₁₀(multiples of free Cmax).** Tested concentrations span up
to 3500-fold. A model linear in raw concentration is dominated by the top point
and its slope is not interpretable — fitted that way, quinidine returns a
*negative* slope alongside a +129 ms effect at 1× Cmax, because the response
saturates and then falls as wells drop out. On the log axis the slope is
milliseconds per decade and the intercept is the effect at 1× free Cmax.

**The source data arrive pre-corrected.** `ddFPDc` in the CiPA myocyte workbook
is already rate-, baseline- and vehicle-corrected by the source study. The
package attaches those three corrections as provenance records marked "applied by
the source study" so the modelling layer's `require_steps` check still passes,
and does not re-apply them.

### 3.3 Informative missingness

A well in arrhythmia has no measurable field potential duration. Of 2067 missing
`ddFPDc` values, 1887 are exactly the wells flagged `EAD = 1`, and missingness
rises 0.9 % → 4.7 % → 26.8 % → 60.8 % across concentration indices 1 to 4. The
endpoint therefore disappears precisely where the drug effect is largest, and a
complete-case slope is biased toward zero.

This is not fixable with the published data — the underlying durations do not
exist — so it is **bounded** instead. Every in-vitro fit is reported twice: over
all concentrations, and restricted to indices 1–3 where missingness is far lower.
For the high-risk drugs the two agree closely (dofetilide 128.7 vs 133.1 ms at
1× free Cmax; sotalol 123.2 vs 132.8), which is the evidence that the ranking is
not an artefact of the dropout.

## 3.4 Validating the in-vitro rate correction

The in-vitro cascade cannot be validated end to end on the CiPA myocyte
workbook, because that workbook publishes the endpoint already corrected. The
CiPA MEA pilot supplies the missing check. Supplemental Table 1 of Millard D et
al., *Toxicol Sci* 2018;164:550-562, doi:10.1093/toxsci/kfy110, reports for each
of its 18 participating studies the mean baseline beat period, the mean field
potential duration, and the mean rate-corrected field potential duration. Given
any two, the third follows from the correction formula.

The table is published as an Enhanced Metafile pasted into the supplementary
Word document. `dish2ecg.io.emf_tables` parses the EMR_EXTTEXTOUTW records and
their reference points, so the cell values and the table layout are recovered
exactly; this is reading the file format, not optical character recognition.

Applying `power_correction(FPD, BP, 1/3)`:

| | |
|---|---|
| studies within 5 ms | 17 of 18 |
| median absolute difference | 0.6 ms |
| worst relative error (excluding study 9) | 1.28 % |
| exponent implied by the published values | 0.3244 (Fridericia 0.3333) |
| Spearman(beat-period CV, abs error) | ρ = 0.592, p = 0.010 |

**What the residual is.** The published FPDc is the mean of per-well corrected
values; what a summary table permits computing is the correction applied to the
mean beat period and mean duration. Those are not equal — the correction is
non-linear, so by Jensen's inequality the gap grows with the dispersion of the
beat period. The measured rank correlation between the beat-period coefficient
of variation and the absolute error is 0.592 (p = 0.010), which is that
prediction confirmed. It is an artefact of the summary, not of either
implementation.

**What the residual is not.** Study 9 (CDI cells on the ACEA platform) is 22.0 ms
off with a beat-period coefficient of variation of only 0.061 — the fifth lowest
in the table. Averaging order cannot account for it. We have no explanation; the
discrepancy is asserted in `tests/test_invitro_validation.py` so that it cannot
be silently absorbed by a future change.

**What this does and does not establish.** It establishes that the rate-
correction step of the in-vitro path reproduces an independent published
implementation. It does not validate the baseline and vehicle subtraction steps,
for which no public uncorrected dataset exists, and it says nothing about
whether FPD translates to QT.

## 3.5 Validating the clinical cascade per record

The Vicente 2016 comparison in section 1 is arm-level. ECGCIPA permits a stronger
check, because its ADaM tables publish the study's own derived values next to the
raw intervals. Over its records:

| our value | their value | n | median diff | SD | max abs |
|---|---|---|---|---|---|
| Fridericia QTcF | `QTCF` | 5740 | 0.000 ms | 0.190 ms | 14.27 ms |
| `jtpeakc_power` | `JTPC` | 5740 | ~0 by construction | 0.15 ms | 11.26 ms |
| baseline subtraction | `CHG` | 1915 | 0.043 ms | 0.218 ms | 5.97 ms |

Only 3 of 5740 QTcF records differ by more than 0.05 ms. The rate correction and
the baseline step therefore reproduce a regulatory-grade derivation essentially
exactly, which is a stronger statement than the arm-level agreement in section 1.

One structural note. ADEG carries the three replicate ECGs per timepoint *and*
the study's own average of them (`EGREPNUM` 0 after filling). Taking that row
rather than re-averaging keeps the analysis identical to the study's. Their
average-row QTcF is the mean of the per-replicate QTcF, whereas Fridericia of the
averaged QT and RR differs from it by about 0.26 ms — Jensen's inequality again,
and small enough to ignore here, but it is the reason the two are not bit-identical.

## 3.6 How much the clinical anchor moves between studies

Ranolazine, verapamil and dofetilide appear in both ECGRDVQ (22-subject crossover,
each subject's own placebo period subtracted) and ECGCIPA. ECGCIPA's design is
mixed (Vicente 2019): ranolazine and verapamil are parallel arms of 10 subjects each;
dofetilide is a part 2 crossover in a further 10 subjects, with no placebo period of
its own. One part 1 placebo group (10 subjects) is subtracted, time-matched, from
**every** ECGCIPA arm, so placebo subtraction is unpaired for all three drugs.

Fitting each study separately and evaluating both at a *matched* plasma
concentration — the lower of the two geometric-mean Cmax values, so neither model
extrapolates:

| drug | matched conc | ECGRDVQ (paired placebo) | ECGCIPA (unpaired placebo) | difference |
|---|---|---|---|---|
| Ranolazine | 2043 ng/mL | 11.2 [7.6, 14.8] | 19.4 [13.1, 25.8] | +8.2 ms |
| Verapamil | 114 ng/mL | 3.5 [-0.1, 7.0] | 8.4 [6.8, 10.0] | +5.0 ms |
| Dofetilide | 1.32 ng/mL | 34.1 [29.8, 38.3] | 43.7 [39.0, 48.5] | +9.7 ms |

ECGCIPA reads higher every time, by 5.0 to 9.7 ms, and for dofetilide the 90 %
intervals do not overlap. **A unit conversion matters here too** — ECGCIPA reports
dofetilide in pg/mL and every other analyte in ng/mL, so comparing across studies
without harmonising silently mixes scales by a factor of 1000.

**How much of the shift is the placebo group, measured.**
`scripts/placebo_commonmode_test.py` refits both studies with *no* placebo
subtraction (baseline-corrected change only) and reports the placebo correction
actually applied:

| drug | ΔΔQTcF shift | shift with no placebo subtraction |
|---|---|---|
| Ranolazine | +8.2 ms | +6.0 ms |
| Verapamil | +5.0 ms | +0.6 ms |
| Dofetilide | +9.7 ms | +7.0 ms |

The mean post-dose placebo change from baseline is −7.6 ms in ECGRDVQ and −9.05 ms in
ECGCIPA's single 10-subject group (−10.1 ms over the records that have a matched
concentration). The shared placebo group therefore accounts for roughly 2–4 ms of
the shift, and for essentially all of it for verapamil. For ranolazine and
dofetilide, **about 6–7 ms remains after placebo subtraction is removed altogether,
and is not explained here.**

An earlier version of this section, and of the README, attributed the whole shift
to unpaired placebo subtraction in a parallel design. That was withdrawn: dofetilide
is not a parallel arm, unpaired subtraction adds variance and shared error without
being expected to bias the estimate, and the measurement above leaves most of the
ranolazine and dofetilide shift standing. Untested candidates for the remainder are
differences in population, ECG extraction and time sampling between the studies,
and pooling of study days in ECGCIPA's multiple-dose arms.

The practical consequence stands: a dish-to-clinic calibration fitted on a single
clinical study inherits an unmodelled between-study error of 5–10 ms for these three
drugs, the same order as the 10 ms regulatory threshold, and no single study can
reveal it.

## 3.7 The ECGDMMLD dofetilide gap, resolved

Earlier versions of this document carried an unreconciled 6-7 ms gap: our
dofetilide-alone ΔΔQTcF of 39.0 ms against a published 45.2 ms. The published value
came from Vicente 2016 Table 3. The **primary** paper for that study, Johannesen L
et al., *Clin Pharmacol Ther* 2015;99:214-223, doi:10.1002/cpt.205, reports its own
numbers in the text:

| quantity | Johannesen 2015 | dish2ecg |
|---|---|---|
| dofetilide alone, ΔΔQTc after the **evening** dose | 39.1 ms | **39.0 ms** |
| dofetilide alone, ΔΔQTc after the **afternoon** dose | 22.4 ms | **22.2 ms** |

Two independent timepoints, both matching to 0.2 ms. The gap was a disagreement
between two published analyses of the same data, not a defect in this pipeline.
Our at-Tmax value is also the maximum over all timepoints (39.0 ms at 13 h), so
45.2 ms is not a maximum-over-time definition either; we cannot account for it.

**A specification difference was tested and rejected as the cause.** Johannesen 2015's
statistical methods specify, for the concentration analysis, "concentration as a fixed
effect and subject as a random effect on intercept **and concentration**" — a random
slope, where this package defaults to a random intercept. `scripts/random_slope_test.py`
refits both crossover studies both ways. The random slope does **not** close the
dofetilide gap, it widens it (-7.7 ms to -11.8 ms against the Vicente 2016 target),
while improving quinidine in ECGRDVQ (-4.0 ms to -1.0 ms). Mean absolute gap over the
nine published comparisons moves from 3.16 ms to 3.04 ms — no material difference. The
intercept-only default is retained and the alternative is available via
`random_slope=True`.

## 3.8 What the waveforms settle: hysteresis, and the delineation rule

ECGRDVQ and ECGDMMLD ship **waveforms** as well as interval tables - `medians/`
(one median beat per record, 16 leads, 1000 Hz, with a `.atr` carrying the study's
own fiducial points) and `raw/` (the 12-lead, 10-second strip the median was built
from). `EGREFID` joins them to the clinical CSV. Earlier versions of this document
stated that only summary intervals were public; that was wrong.

Decoding the `.atr` of a median beat and differencing the fiducials returns the
clinical table's own intervals exactly - over 242 beats the median difference
between the annotation-derived QT and the table QT is **0.00 ms**. That identity is
what licenses everything below.

**Hysteresis correction is not merely unimplemented here - it is unidentifiable.**
QT needs one to two minutes to settle after a rate change, while a record is 10
seconds, so the window covers only **8-17% of one adaptation time constant**: the
RR history that determines the QT in the strip mostly precedes the strip. No
weighting scheme over 10 seconds of RR can recover it.

What that costs is measurable. On 238 records (22 subjects, 5 arms, RR 661-1405 ms),
with R peaks detected from lead II - the detector's mean RR matches the published RR
to a median of **+0.10 ms**, SD 5.6 ms, with |error| > 20 ms in 0.8% of records:

| | |
|---|---|
| beats per 10 s window | median 11 (8-15) |
| within-window RR standard deviation | median 29.4 ms |
| within-window RR range | median 87.0 ms |
| QTcF spread across the window's slowest and fastest beat | median **12.6 ms** |
| records where that spread exceeds 10 ms | **66.8%** |

The spread is an upper bound on the pairing choice, not an error estimate - the
table pairs QT with the window's *mean* RR, which sits in the middle of that range,
and is the right choice. The point is that the quantity being averaged moves by more
than the regulatory threshold inside a single recording, in two thirds of records.
This is sinus arrhythmia in young healthy adults, not noise.

**The T-end rule is worth a few milliseconds, and the difference is rate-dependent.**
An independent tangent (Lepeschkin) T-end measured on the same median beats, against
the study's own annotation (n = 220 after rejecting 9.1% of tangent fits as
implausible - a limitation of this simple implementation, not of the data):

| | |
|---|---|
| tangent minus study T-end | median **-4.6 ms**, IQR [-7.3, -3.0] |
| rate dependence | **-0.0317 ms per ms RR**, i.e. -9.5 ms across a 300 ms RR change (p = 0.026) |

A *constant* offset between delineation rules cancels in the ΔΔ cascade. A
**rate-dependent** one does not, because heart rate moves between baseline and
on-drug. This reproduces, on this corpus and with a different algorithm pair, the
direction and order of magnitude found on PTB-XL+ (-5.9 ms across 300 ms RR between
two commercial algorithms). Two independent datasets, two different comparisons, the
same sign and scale - which makes extraction method a real effect on *absolute*
QT comparisons across algorithms.

It is **not** an explanation of the between-study shift of section 3.6. That was
tested directly (docs/DELINEATION_AND_SHIFT.md): because the cascade is linear in
QT, the movement it produces equals the cascade applied to the per-record
difference. On 3,163 ECGRDVQ and 2,515 ECGCIPA median beats the baseline- and
placebo-corrected movement is -1.8 to +0.3 ms, and the *between-study* difference
in it is -0.32 ms (ranolazine) and -2.08 ms (dofetilide) against a 6-7 ms
discrepancy. The reason is exposure: the placebo-corrected within-subject RR
change in these trials is only -33 to +18 ms, not the 300 ms span used for the
sensitivity above, so the predicted effect is ~1 ms and the observed one matches.
The 6-7 ms shift remains unexplained.

## 4. The bridge

No conversion of FPD into QT is attempted or permitted. Comparison is on three
shared scales:

1. **Concentration** — multiples of the free (unbound) clinical Cmax. Clinical
   plasma concentrations are converted ng/mL → mol/L using the PubChem molecular
   weight, multiplied by the fraction unbound, and divided by the free Cmax used
   by the in-vitro study.

   **Source of the free Cmax, and why the choice does not matter.** Crumb WJ et al.,
   *J Pharmacol Toxicol Methods* 2016;81:251-262, doi:10.1016/j.vascn.2016.03.009,
   Table 1 tabulates free Cmax in nM for 30 drugs; it is reproduced in
   `data/crumb2016_free_cmax.csv`. For the 12 drugs shared with the CiPA myocyte
   workbook, Crumb and the Blinova reference agree to a median ratio of 1.000 and
   11 of 12 fall within 15 % — quinidine is the sole exception (Crumb 842.9 nM
   against Blinova 3000 nM). The CiPA in-silico training table diverges further
   (terfenadine 13.3x, quinidine 3.8x, verapamil 1.8x). **The bridge is invariant to
   this choice**: the same reference divides both the clinical and the in-vitro
   concentration, so a change of reference rescales both axes identically and leaves
   every comparison, and the `point_inside_dish_range` flag, unchanged. The reference
   does matter for any statement of the form "percent channel block at 1x free Cmax",
   which is why `results/channel_block_at_crumb_free_cmax.csv` recomputes those on
   Crumb's own values — internally consistent with the IC50s measured in the same paper.
2. **Effect** — percent change from the corrected baseline, or a standardized
   effect. Absolute milliseconds are reported side by side but are explicitly not
   claimed to be commensurable.
3. **Decision** — whether the upper 90 % bound crosses a threshold. The clinical
   threshold is the 10 ms of ICH E14. There is **no** agreed in-vitro equivalent,
   so `decision_agreement` requires the caller to supply one and records that it
   is a user convention, not a regulatory criterion.

**The comparison point.** Both sides are evaluated at the geometric-mean
*observed* clinical free Cmax, not at a nominal 1×. In ECGRDVQ the observed
exposure reaches only 0.29× the free Cmax for quinidine and 0.50× for verapamil;
evaluating the clinical model at 1× would extrapolate 3.5-fold beyond the
measured range. The quinidine comparison point of 0.29× sits just below the
in-vitro tested range of 0.32–10×, which the bridge table flags with
`point_inside_dish_range = False`.

## 4.1 What the bridge is not

The bridge relates a drug's in-vitro concentration-response to that drug's
clinical concentration-QTc relationship. It is **not** a personalised predictor,
and the distinction is empirical rather than cautionary. Blinova 2019
(doi:10.1111/cts.12674) ran the subject-matched experiment: iPSC-CMs derived from
16 of the same healthy volunteers who supplied clinical concentration-QTc data,
dosed with dofetilide and moxifloxacin. Subject-specific ΔΔAPD₉₀c slopes showed
no rank correlation with the same subjects' clinical ΔΔQTc slopes (ρ = −0.106,
p = 0.70 for dofetilide; ρ = −0.089, p = 0.75 for moxifloxacin), and the
null survived switching to Fridericia-corrected ΔΔAPD₉₀cF.

The informative control in that paper is that the subjects' clinical ΔΔQTc
responses to the two drugs did not correlate with each other either (ρ = 0.096,
p = 0.73). Between-subject variation in QT response is therefore not a stable
personal trait at this sample size, so the in-vitro null is not by itself
evidence that the dish is the weak link.

## 5. Data provenance

Every derived value carries an ordered tuple of `CorrectionRecord`s — step,
method, the formula actually applied, the parameters, the citation, and any note
(such as how many records fell back to a pooled exponent). `Endpoint.describe()`
prints the cascade; `provenance_columns()` flattens it for a results table.
Slicing an endpoint uses `Endpoint.subset(mask)`, which keeps the cascade; slicing
`.values` directly would strip it.

## 6. Known data defects repaired

- **ECGRDVQ PR wrap.** Two PR values are stored as 32-bit wrapped negatives
  (−4294966951 and −4294966972; subject 1007, Verapamil HCL, TPT 1.0). Adding 2³²
  recovers 345 and 324 ms, consistent with first-degree AV block under verapamil.
  An unguarded mean of PR is −1 644 475 ms.
- **TPEAKTPEAKP** is populated for 42 of 5232 ECGRDVQ records and none of
  ECGDMMLD. Unusable; not analysed.
- **Concentration units are not uniform.** Dofetilide is reported in pg/mL in
  both studies, every other analyte in ng/mL. Harmonised at load.
- **Analyte and treatment names differ** (`Quinidine` vs `Quinidine Sulph`,
  `Verapamil` vs `Verapamil HCL`). Joined on a normalised key.
- **Drug-name spelling in the CiPA workbook.** `D,l Sotalol` and `D,l,Sotalol`
  are one compound under two spellings; normalised, and the normalisation is
  flagged per row.
- **Micro sign in the supplement tables.** The per-drug concentration tables are
  headed "Drug Concentration, M" but the values are micromolar; the micro sign is
  a Symbol-font glyph that extracts as a private-use codepoint. All molar values
  in the lookup are the printed value × 1e−6.

## 7. Not assessed

Residual diagnostics, normality of the random effects, and influence analysis
were not performed. Hysteresis correction of RR was not applied, for want of
beat-to-beat RR history in the public tables. The nine ECGRDVQ raw records with no
median-beat counterpart were not investigated; waveform files were not used at
all.
