# Changelog

## 0.2.2 — 2026-11-04

Fraction unbound for bepridil, chlorpromazine, cisapride and terfenadine
(`data/drug_reference.csv`, `docs/DRUG_REFERENCE_SOURCES.md`).

- Values: bepridil 0.005 [0.003-0.010], chlorpromazine 0.055 [0.03-0.08], cisapride
  0.0225 [0.020-0.025], terfenadine 0.03. From product labels (reproduced on
  drug-information sites; none of the four has a current SPL) and, for bepridil, a
  measured equilibrium-dialysis study (Pritchard 1985, >99.7% bound). Bepridil's value is a
  bounded estimate: its binding is saturable, so fu is concentration-dependent.
- Cross-checks: bepridil fu implied by Crumb free Cmax over the label's total Cmax is 0.0050;
  chlorpromazine's implied total sits inside the range the label binding was measured over;
  cisapride is in mild tension (published peaks imply 0.009-0.017, label gives 0.020-0.025).
- **Correction to an earlier statement in this repository**: the provenance note said Crumb 2016
  and Li 2017 tabulate protein binding for this drug set. Neither does (checked in the PDFs).
- **Terfenadine's CiPA free Cmax (4 nM) is not comparable to Crumb/Blinova's (0.3 nM)**: with fu
  0.03 they imply ~63 vs ~4.7 ng/mL total, against a reported normal-dose mean peak of ~1.5
  ng/mL, and a Redfern 2003 excerpt says 4 nM was reached with a P450 inhibitor present. This is
  the likely cause of the 13x divergence between the two tables.
- None of these four drugs is in the dish-clinic bridge; no published result changes.
- 7 tests added.
- Added `docs/PROJECT_REPORT.md`, a narrative research report covering the whole project (also issued as a
  self-contained HTML edition with figures embedded).
- Added `results/pipeline_diagram.png` (`scripts/make_pipeline_diagram.py`) as Figure 2: the in-vitro, in-vivo and
  in-silico pipelines against the shared cascade steps, with the FPD-to-QT guardrail drawn where it fires. Every
  panel below the top row plots this repository's own data (18 Millard studies FPD->FPDc; 439 drug-free clinical
  records; the dofetilide baseline/placebo waterfall at 2.5 h; the fitted in-vitro line; 328 clinical records; the
  ORd dose-response). Later figures renumbered 3-6.
- **Found while drawing the figure:** in `results/ord_apd90.csv` the `ok` flag is true for all 408 runs, but beyond
  about 3.8x free Cmax the model stops repolarising - APD90 is NaN, or is measured from a depolarised resting
  potential (v_rest > -60 mV) and is meaningless. The flag does not capture repolarisation failure; filter on the
  physiology. The affected rows were never used in a reported result.
- Report sections 10 (challenges met and what they cost) and 11 (questions for the study authors) added; Availability
  renumbered to 12. Section 11 is also shipped standalone as `docs/QUESTIONS_FOR_AUTHORS.md`.
- `scripts/relative_scale_test.py` -> `results/relative_scale_test.csv`. Tests whether normalising removes the
  dish-clinic gap. On a percent-of-baseline scale the gap narrows from 1.81x to 1.48x but does not close; relative to
  a reference drug the two PURE hERG blockers agree almost exactly (dofetilide-normalised quinidine: dish 1.023 vs
  clinic 1.025, ratio 0.999), while ranolazine is 1.92x off and verapamil inverts. The dish baseline needed for a
  percent scale is NOT in the Blinova workbook and is imported from Millard 2018 (median control FPDc 482 ms); this is
  flagged in the output. Verapamil's clinical interval spans zero, so its ratio has unbounded variance.
- `scripts/bland_altman.py` -> `results/bland_altman_stats.csv`, `results/bland_altman.png`. Method comparison where
  it is licensed: optical APD90 against electrode FPD (97 paired drug x concentration cells) and our tangent T-end
  against the study's (220 records, 22 subjects). Optical vs electrode: bias +34.7 ms, limits of agreement -227 to
  +297 ms, proportional-bias slope 0.800 (p < 0.001; still 0.502 with the one extreme cell removed) - the two are not
  interchangeable, which is the sharpest available test of the source paper's "equivalent measures" wording. Tangent
  vs study: bias -9.1 ms, limits -64 to +46 ms, no proportional bias (p = 0.49).
  Limits of agreement use the closed-form unbalanced one-way variance components (Bland & Altman's repeated-measures
  treatment) because the mixed-model fit did not converge on the clustered optical comparison; the naive SD of
  differences would have understated the spread.
- Figure 2 extended: 90 % confidence bands on all three concentration-response panels (dish and clinic from the
  fitted SEs with a t quantile on n_groups-1, matching the cascade's own intervals; in silico as the range over the
  three cell types), the ICH E14 10 ms line on the clinical panel and on the comparison panel, and the per-drug E14
  call-out (dofetilide 77.7, quinidine 80.5, ranolazine 14.8 above; verapamil 7.0 below). The dish panel states that
  no in-vitro equivalent of the 10 ms rule exists.
- Report sections 12-17 added: the biology of the blockers and how the three signals arise; the four corrections
  (rate, hysteresis, autonomic tone, irregular beats and the median beat) with what each is doing and why hysteresis
  is unidentifiable here; why one solution does not fit all three arms; whether cross-system consistency is a
  sensible question, with the diagnostic-accuracy and Bayesian framing from the translation literature; what the data
  cannot support (plate-level and non-hERG channel random effects) and a ranked list of the data we would want; and a
  per-paper literature summary with the three papers read only via abstract flagged as such. Availability renumbered
  to 18; 13 references added.
- `dish2ecg.viz` gains the uncertainty layer, so bands and agreement limits are library code rather than
  script-local: `confidence_band`, `plot_confidence_band`, `variance_components`, `bland_altman`,
  `plot_bland_altman`. The band reproduces the interval the model layer reports (t on n_groups-1), asserted in a test
  against the shipped in-vitro interval, so a plotted band and a quoted effect cannot disagree. 9 tests added.
- `scripts/confidence_bands.py` -> `results/confidence_bands.png`, `results/confidence_bands.csv`: one panel per arm
  for dofetilide, each read at the geometric-mean observed free Cmax and annotated with the interval there
  (in vitro +132 [93, 171]; in vivo +73 [69, 78] with the 10 ms line; in silico +50 [41, 50], labelled as a range
  over cell types and NOT a confidence interval).
- `scripts/methodology_sensitivity.py` -> `results/methodology_sensitivity.csv`: every processing choice refitted
  against its alternative. **Bazett flips the ICH E14 decision for verapamil** (Fridericia 3.5 [-0.1, 7.0], negative;
  Bazett 10.5 [6.0, 15.0], positive), now asserted in a test. Forcing a linear in-vitro predictor halves verapamil
  (-87 to -43) and cuts dofetilide by a quarter (+132 to +96). Dataset as a fixed rather than random effect moves
  dofetilide 13.7 ms. The white-paper concentration-QTc model differs from our pre-formed double delta by at most
  3.4 ms with no decision changes; the missing centred-baseline covariate is named as the deviation to fix first.
- Report section 18 (methodological audit) added, with each mixed model written out; Availability renumbered to 19.
- `scripts/model_shape.py` -> `results/model_shape.csv`: linear vs log-linear vs Emax on the clinical arm, same
  random intercept, ML so BIC is comparable, shape parameters profiled. Linear is best for dofetilide only;
  quinidine prefers log-linear (dBIC 19.3) and verapamil Emax (dBIC 7.1). The effect at Cmax moves at most 6 ms and
  NO ICH E14 decision changes - unlike the rate-correction choice, which flips one.
- Report section 18.7: the in-silico arm fits no concentration-response model at all, and the curvature in
  Figure 2c is largely a log-axis artefact (linear fit of dAPD90 on multiples gives R2 = 0.9993 dofetilide,
  0.9894 quinidine, 0.9820 ranolazine, 0.9089 verapamil; dofetilide endo extended to 30x falls to 0.9387).
- Repository URL https://github.com/rejaul/dish2ecg recorded in README, pyproject [project.urls] and CITATION.cff;
  the archive now carries its own git history with that remote configured.
- Licence and authorship set: MIT, Copyright (c) 2026 Rejaul Karim. `CITATION.cff` gains the top-level `authors`
  block that CFF 1.2.0 requires and that was missing, and `pyproject.toml` gains an `authors` entry.
- **Citation fix in `CITATION.cff`:** the Scientific Reports 2019 J-to-Tpeak rate-correction paper
  (doi:10.1038/s41598-019-51491-4) was attributed to Vicente; the first author is Hnatkova. Corrected.
- `CITATION.cff` version bumped from 0.1.0 to 0.2.2 and the release date to 2026-11-06.
- `docs/PROJECT_REPORT.md` rewritten in plain prose with the same 19 sections, and Figures 7 (Bland-Altman) and 8
  (confidence bands) added. Every number was compared against the previous version; each one added or dropped is
  accounted for (figure-8 values, Park 2018's reported sensitivities, the 142 test count, journal page numbers).
- **Errors found while rewriting, and fixed:**
  - The 408 ORd simulations were described as "the four bridge drugs and 34 concentrations". They cover the 12
    CiPA training drugs, three cell types and 11 or 12 concentrations each; four of the 12 are bridge drugs. The
    Figure 2 caption, which merged the 12-drug APD90 set and the 28-drug qNet set, now names both.
  - Bayesian post-test probabilities (20 % pre-test, 53 % / 9 % / 4 % after hERG results) were quoted in the report
    and in a chat summary and credited to Gintant 2022. They are not in Gintant 2022 (a single-author review of
    diagnostic-performance principles; abstract only), nor in Strauss 2021 or Park 2018, both read in full.
    Strauss 2021 attributes the 20 % pre-test figure to Vargas et al. and says it disagrees with their exact
    probabilities. The figures are removed; what was verified replaces them.
  - Park 2018 (first author Park, with Gintant as a co-author) is the source of the sensitivity and specificity
    figures. Its action potential assay is the conventional non-clinical one, not hiPSC-CM.
  - Of the 136 drug-by-concentration conditions, the M cell is the longest APD90 in 121. All 15 exceptions are
    conditions where the M cell itself failed to repolarise, so it is the longest in every condition where it
    repolarises normally.
  - The reference list no longer carries Vicente 2018 (cited nowhere in the text). Garnett 2018 is added, and four
    sources seen only as abstracts or search excerpts are marked as such.
- **Doc fix:** the bootstrap intervals in the held-out risk table (`docs/SUMMARY.md`, `docs/RISK_CLASSIFICATION_28.md`)
  were stale (for example qNet Spearman interval [0.55, 0.95]). `scripts/risk_comparison.py` was re-run from the wells
  table and reproduces `results/risk28_metrics_ci.csv` exactly; the docs now quote that file (qNet [0.52, 0.95]).
  Point estimates never changed. The all-28 qNet AUC-low interval was corrected from [0.86, 1.00] to [0.85, 1.00].
- `docs/RESEARCH_REPORT.md` carries a banner: frozen at 0.1.1, listing the four statements that have since changed.
- `docs/PROJECT_REPORT.md`: the limits section now treats the dish exposure axis properly. Adsorption to plate
  plastic and tubing, and partition into the cell membrane and intracellular lipid, are named as a loss distinct from
  serum protein binding, with quantified cardiac examples verified against the source full text (Alvarez Baron 2025,
  doi:10.1111/cts.70451: LC-MS/MS concentration verification in hERG rigs gave losses of 6% ranolazine, 19%
  chloroquine, 41% verapamil, 58% vanoxerine; correcting moved the hERG IC50 to 2.2 nM for vanoxerine and 0.75 uM for
  chloroquine) and the direction of that one mechanism (free < nominal shifts dish curves toward higher potency, so the
  over-read would widen). Two of those compounds are bridge drugs and their losses differ sevenfold, so the loss does
  not cancel when drugs are compared on a shared exposure axis.
  **Correction within this release:** the first draft of this paragraph stated chloroquine lost "about 2%" in the same
  system. That conflated the ICa/INa UV-spectroscopy figures (2.5% and 2.2%) with the hERG LC-MS/MS figure (19%), and
  was written before the source full text was read. Corrected after verification. Explicitly NOT a reinstatement of withdrawn claim 1: the net sign of a
  free-concentration correction is still unestablished (Schocken 2018 runs the other way). A measurement of free
  compound in the well was added to the "what would settle it" list, plus a glossary entry for nominal concentration.
- **Chlorpromazine binding citation read at abstract level**: the note that chlorpromazine's free fraction falls as alpha-1-acid
  glycoprotein rises had rested on a truncated excerpt. The PubMed abstract (Verbeeck, Cardinal and Wallace, Eur J Clin Pharmacol
  1984;27:91-97, PMID 6489431) confirms it: 64 healthy volunteers, negative correlation, p<0.005. `data/drug_reference.csv` now names
  the paper. The full text is closed, so the value of fu itself (0.055) still comes from the product label, as before.
- **Correction to the quinidine/hysteresis citation in `docs/PROJECT_REPORT.md`**: the report quoted a BIC improvement of
  "about 110" from hysteresis and a "gap of about 4" between Emax and linear, attributed to an unnamed paper. A later
  search showed the paper (Cellière, Krause, Bonnefois, Chauvin; PMID 40347307, PMC12065759). Its hysteresis table gives
  105.8 for ΔΔQTc and 109.6 for ΔQTc. The Emax-versus-linear gap is not visible in the excerpt we read, so that number is
  removed, and the conclusion drawn from comparing the two is now stated as our own judgement. The reference now names the
  authors and sits in alphabetical order. The paper is still known only from a search excerpt. `docs/RESEARCH_REPORT.md`
  is frozen and still says the authors were not retrieved.


## 0.2.1 — 2026-11-04

Reproducibility audit.

- **14 of the 38 result files had no producing script** — they came from ad-hoc
  analysis and could not be regenerated. Of those: **9 are now reproducible**
  (`scripts/risk_comparison.py` and `scripts/delineation_shift.py` added, both
  verified to reproduce the committed CSVs exactly at max difference 0.00, and
  `arm_confounding.py` extended to write `arm_perdrug.csv`); **1 was deleted**
  (`results/waveform_qtc_window.csv`, a redundant superset of the
  script-produced `waveform_rr_window.csv`); and **4 remain ad-hoc** and are
  labelled as such in `results/README.md` with their inputs named. The other 24
  files already had a producer.
- `results/README.md` added: every result file mapped to the script that writes
  it, with those four ad-hoc ones labelled and their inputs named, plus
  the full regeneration order.
- Test count claims corrected (117 -> 126).

## 0.2.0 — 2026-11-04

Two results that each overturn an earlier claim of this repository.

- **Out-of-sample TdP risk test, all 28 CiPA drugs** (docs/RISK_CLASSIFICATION_28.md).
  The validation-set parameters are public on the FDA/CiPA `Model-Validation-2018`
  branch; those 28 drugs are exactly the Blinova panel. 2,240 simulations, 0
  failures, hERG uncertainty carried by a 20-sample bootstrap ensemble per drug
  and summarised in output space (validated against the published point fit on
  the 12 training drugs: within 0.0006 uC/uF for 11 of 12).
  qNet separates high risk perfectly on all 28 including the 16 held out
  (AUC 1.000). **WITHDRAWN**: the claim that the dish "detects danger but does not
  grade safety" was an artifact of the 12-drug training subset - on the held-out
  16 the dish's AUC for low risk is 0.945, against qNet's 0.909.
- **The T-end rule does not explain the between-study shift**
  (docs/DELINEATION_AND_SHIFT.md). 3,163 ECGRDVQ + 2,515 ECGCIPA median beats
  delineated with identical code. Per record the rules differ by a median -4.5
  and -7.1 ms, but the baseline- and placebo-corrected movement is -1.8 to
  +0.3 ms and the between-study difference is -0.32 ms (ranolazine) / -2.08 ms
  (dofetilide) against 6-7 ms. **WITHDRAWN**: "extraction method is a live
  contributor to the unexplained between-study shift" was too strong. The
  rate-dependence is real but the placebo-corrected RR exposure in these trials
  is only -33 to +18 ms, not 300 ms, so it is worth ~1 ms. The shift remains
  unexplained.
- `viz.read_wfdb` extended to WFDB format 32 and to gains quoted per uV
  (ECGCIPA medians); signals are now always returned in mV. It still refuses
  unknown formats rather than guessing.
- 9 tests added.

## 0.1.9 — 2026-11-04

- Corrected an overstatement attributed to Blinova 2018. The paper names **three**
  most-useful predictors (arrhythmia-like events at any concentration, prolongation
  at any concentration, and prolongation at clinical Cmax); it does not single out
  event counts as "the best". Fixed in SUMMARY.md and RISK_CLASSIFICATION.md.
- Added an explicit source block to SUMMARY.md for the two phrases quoted from
  Blinova K et al., Cell Reports 2018;24(13):3582-3592,
  doi:10.1016/j.celrep.2018.08.079 — the quotes are from the paper, not from the
  published well-level data file that carries the same name.

## 0.1.8 — 2026-11-04

- `docs/SUMMARY.md`: a reviewer-facing account of the whole project — what the
  bridge shows, what the assay is good at, what cannot be estimated and why, the
  validation record, and an explicit contribution statement.
- Removed `results/waveform_delineation_ok.csv`, a stale orphan from the ad-hoc
  waveform run (215 rows, median -4.55 ms) that contradicted the script's own
  output (`waveform_delineation.csv`, 220 in-range of 242, median -4.63 ms) and
  was referenced by nothing.
- Corrected the Millard 2018 rate-correction agreement from "median 0.6 ms" to
  **0.65 ms** in four documents; the value read back from
  `results/invitro_correction_validation.csv` is 0.65.

## 0.1.7 — 2026-11-04

> Its headline conclusion about the dish was reversed in 0.2.0 by the 28-drug
> out-of-sample test; read that entry alongside this one.

CiPA qNet and a TdP-risk classification comparison (docs/RISK_CLASSIFICATION.md).

- `scripts/qnet_simulation.py`: ORd-CiPA v1.0 (Dutta 2017) with the model's own
  dynamic hERG Markov binding plus Hill block on the other six channels, at the
  CiPA protocol (endocardial, 0.5 Hz, 1-4x therapeutic). 48 simulations, 0
  failures. Control qNet 0.0785 uC/uF and APD90 286.7 ms, both matching published
  ORd-CiPA values; dofetilide at 4x prolongs APD90 to 416.8 ms.
- qNet puts the three CiPA risk categories in completely non-overlapping bands
  (low 0.0831-0.1009, intermediate 0.0687-0.0747, high 0.0390-0.0648),
  reproducing the published training-set separation.
- **This is IN-SAMPLE**: those 12 drugs are the set qNet was optimised on, so the
  comparison does not show the simulation beating the dish. The out-of-sample
  finding is the shape of the dish's competence: its arrhythmia-event fraction
  identifies all four high-risk drugs (AUC 1.00) but separates low from the rest
  at only AUC 0.69 [0.36, 1.00]. The dish detects danger; it does not grade safety.
- Three traps recorded: qNet excludes INa and uses only ICaL's calcium carrier
  (including INa gives 0.0044, an 18x error that still looks plausible); Kt is a
  fixed 3.5e-5 that the mmt defaults to 0, which halves dofetilide's effect; and
  ikr.D must carry the nanomolar number despite its [mM] annotation.
- 6 tests added.

## 0.1.6 — 2026-11-04

The APD90 -> FPD -> QT ladder (docs/APD90_FPD_QT_LADDER.md).

- **Arrow 1 (APD90 -> FPD) is not identifiable in public data.** Blinova 2018
  pools both endpoints into one column; the only optical (APD90) platform is
  CellOPTIQ at site 2, which is also the only serum-free site. Leave-one-site-out
  over 112 drug x concentration cells: the optical arm's gain vs the other nine
  sites is 1.436 +- 0.146, the highest of all ten (z = +2.16 vs the MEA sites,
  CI excluding 1.0); paired, optical = 1.51 x MEA + 24.8 ms (n = 97). Endpoint,
  site and medium cannot be separated, and Schocken 2018 shows the medium term
  is both large (104 ms baseline) and of unpredictable sign.
- **In-silico rung added**: O'Hara-Rudy 2011 + CiPA multichannel IC50s,
  408 simulations, 0 failures, controls matching the published APD90s
  (endo 265.5 / epi 228.0 / mid 329.9 ms). At the bridge match points it agrees
  with the clinic on 4/4 signs where the dish manages 3/4 - including verapamil,
  which the dish inverts by 86 ms. Magnitudes are not calibrated (0.60-3.5x).
- **Arrow 3 is a coupling problem, not an average**: with uncoupled cells the
  QT surrogate is just the M-cell APD90 - the M cell repolarises last in 121 of
  136 conditions, including every therapeutic match point. The 15 exceptions are
  all at >= 3.78x free Cmax, where the M cell loses a stable single
  repolarisation (two fail outright) and the surrogate falls back to endo. The coupled 1D cable is written
  (scripts/cable_simulation.py) but blocked on this host - pocl advertises
  cl_khr_fp64 on the device but not the platform, so myokit refuses double
  precision and ORd is unstable in single. Also documents the myokit trap that a
  single-cell model file has no diffusion_current binding, so set_conductance()
  silently does nothing.
- 7 tests added.

## 0.1.5 — 2026-11-04

PhysioNet ECGRDVQ/ECGDMMLD ship WAVEFORMS (`medians/`, `raw/`), joined to the
clinical CSV by `EGREFID`. Earlier docs said only intervals were public - wrong.

- `scripts/waveform_analysis.py` + `results/waveform_*.csv`.
- **Hysteresis is unidentifiable, not just unimplemented**: a record is 10 s,
  which is 8-17% of one QT adaptation time constant. Within-window RR has a median
  SD of 29.4 ms and median range 87.0 ms; the implied QTcF spread across the
  window's extreme beats is a median 12.6 ms and exceeds 10 ms in 66.8% of 238
  records. R-detector validates against the published RR to a median of +0.10 ms.
- **Delineation rule is rate-dependent**: an independent tangent T-end differs from
  the study's by a median -4.6 ms, with a slope of -0.0317 ms per ms RR
  (-9.5 ms over a 300 ms RR change, p = 0.026). Same sign and order as the PTB-XL+
  result, so a constant-offset assumption is wrong and the difference does not
  cancel in the cascade.
- Annotation decoding validated: annotation-derived QT matches the clinical table
  to a median of 0.00 ms over 242 beats.
- 6 tests added. Fixed `ndarray.ptp()` (removed in NumPy 2.0) and replaced silent
  `except: pass` with counted failures.

## 0.1.4 — 2026-11-04

- **`dish2ecg.viz`**: visualisation of what the cascade measures. Minimal WFDB
  format-16 reader and annotation decoder (no new dependency), `plot_ecg_beat`,
  `plot_field_potential`, `plot_qt_rr`, `plot_concentration_response`.
- `scripts/make_figures.py` writes three figures. The ECG panel is a REAL median
  beat from PhysioNet ECGRDVQ, annotated with the study's own fiducial points;
  the decoded fiducial gaps reproduce the clinical table's PR/QRS/J-Tpeak/
  Tpeak-Tend/QT exactly, which the test suite asserts. The field-potential panel
  is labelled a schematic because no public dataset here ships hiPSC-CM waveforms.
- 7 tests added (98 total).

## 0.1.3 — 2026-11-04

- **docs/TRANSLATION_REVIEW.md**: literature review of why dish and clinic differ.
  28 PubMed queries (206 records) plus targeted retrieval, 34 works cited, all
  CrossRef-verified. 20 mechanisms tabulated with direction of bias.
  Notable findings: autonomic tone *attenuates* drug-induced QT prolongation
  (Smith 2007: ibutilide +43 ms with intact tone vs +68 ms additional under double
  autonomic block), so a denervated dish is expected to over-read — a mechanism
  distinct from immaturity and apparently not invoked in the hiPSC-CM translation
  literature. Serum effects contradict a pure protein-binding account (Schocken 2018).
  At least three mechanisms have unknown sign, so the 1.8x inflation factor must not
  be used as a calibration constant.
- REFERENCES.md extended with the review's 30 new works.

## 0.1.2 — 2026-11-04

Seven previously paywalled or unmined papers were supplied and read.

- **The ECGDMMLD dofetilide gap is CLOSED and was never ours.** Johannesen 2015
  (doi:10.1002/cpt.205) reports dofetilide alone at 39.1 ms (evening) and 22.4 ms
  (afternoon); this pipeline gives 39.0 and 22.2. The 6-7 ms gap was against Vicente
  2016 Table 3 (45.2 ms), which disagrees with its own study's primary paper.
- **Random-slope specification tested and rejected.** Johannesen 2015 specifies a
  random effect on intercept *and* concentration. `scripts/random_slope_test.py`:
  mean absolute gap over nine published comparisons 3.16 -> 3.04 ms; dofetilide
  worsens, quinidine improves. Intercept-only retained.
- **J-Tpeak exponent confirmed at source.** Johannesen 2014 (doi:10.1038/clpt.2013.238)
  states the coefficient was derived in that study by pooling 34 thorough-QT studies.
- **Crumb 2016 Table 1 free Cmax recovered** for 30 drugs (`data/crumb2016_free_cmax.csv`).
  Matches the Blinova reference for 11 of 12 shared drugs (median ratio 1.000);
  quinidine is the exception. Channel-block percentages recomputed on it.
  Documented: the bridge is invariant to this choice, since the same reference
  scales both axes.
- **Schocken 2018 mined.** Medium alone moves baseline FPDc 575 ms (serum) vs 471 ms
  (serum-free). 8 of 18 arrhythmia-inducing drugs did so at *lower* concentrations
  with serum, against what protein binding alone predicts; the earlier claim that
  free-concentration correction would simply widen the dish-clinic discrepancy is
  withdrawn.
- 5 tests added.

## 0.1.1 — 2026-11-04

Corrections found while preparing the research report.

- **Cross-study shift mis-explained.** The replication table labelled dofetilide's
  ECGCIPA arm "parallel"; it is a part 2 crossover. The explanation offered
  (unpaired placebo subtraction in a parallel design) is withdrawn. New diagnostic
  `scripts/placebo_commonmode_test.py`: the shared placebo group explains about 2–4 ms
  of the 5.0–9.7 ms shift, all of it for verapamil; 6–7 ms remains for ranolazine
  and dofetilide, unexplained. Replication table regenerated with correct design
  labels; numbers unchanged.
- **J-Tpeak exponent is published.** Vicente 2019 states J-Tpeakc = J-Tpeak/(RR/1000)^0.58
  (attributed to Johannesen 2014). The 0.5799 recovered from PhysioNet is a
  confirmation, not a discovery; citation strings and docs corrected. Swapping the
  power law into ECGRDVQ moves ΔΔJ-Tpeakc by at most 1.2 ms.
- README no longer claims clinical concentration-QTc analysis has "no open
  implementation" (an R package driving MonolixSuite, and R power code, exist) and no
  longer names "PhysioMEA", which could not be verified.
- REFERENCES.md extended; every addition resolved against CrossRef.
- Appendix 8.1: eight public ECG corpora assessed; none carries drug exposure. PTB-XL+
  yielded a measurement of interval-extraction variability (uncertainty 9b):
  rate-dependent QT bias of ~6 ms across a 300 ms RR change between two commercial
  algorithms, and a fitted QT/RR exponent of 0.3355 vs 0.3734 on identical records.
- 5 tests pin the diagnostic and the corrected labels: 86 tests in all
  (`pytest --collect-only -q`).

## 0.1.0 — 2026-11-02

First release.

- Rate corrections: Bazett, Fridericia, Framingham, Hodges, generic power,
  pooled and per-subject fitted exponents, and the Hnatkova 2019 J-Tpeak
  correction.
- Baseline / vehicle / placebo subtraction cascade shared by both domains, with
  enforced ordering and counted unmatched records.
- Mixed-effects concentration–response model per ICH E14 Q&A, with target
  centring, internal predictor scaling, a documented OLS boundary fallback, and
  t-quantile intervals.
- `Endpoint` provenance objects; in-vitro endpoints refuse to be reported as QT
  intervals.
- PhysioNet loaders with checksum verification; MEA export adapter with an
  explicit alias table and no silent inference.
- Reference analysis of ECGRDVQ and ECGDMMLD, validated against Vicente 2016
  Table 3 (published value inside our 95 % CI for 5 of 6 comparisons).
- In-vitro analysis of the CiPA myocyte well-level data with an informative-
  missingness sensitivity analysis, and a four-drug dish-to-clinic bridge.
- In-vitro rate correction validated against the independently published FPDc of
  the CiPA MEA pilot (Millard 2018): 17 of 18 studies within 5 ms, median 0.65 ms,
  implied exponent 0.3244 against Fridericia's 0.3333. Study 9 does not reproduce
  and the discrepancy is asserted rather than hidden.
- `dish2ecg.io.emf_tables`: recovers tabular text from Enhanced Metafile images,
  for supplementary tables published as pasted vector graphics.
- Third clinical corpus: PhysioNet ECGCIPA (CiPA ECG validation study, 60
  subjects, adds chloroquine and lopinavir/ritonavir). `load_ecgcipa` and
  `ecgcipa_wide` handle its CDISC ADaM layout.
- `jtpeakc_power`: the FDA power-law J-Tpeak correction, exponent 0.58, recovered
  as 0.5799 from 5740 records of the study's own derived values. Differs from the
  Hnatkova 2019 linear form by up to 21 ms.
- Clinical cascade validated per record against ECGCIPA's own derived values:
  QTcF median difference 0.000 ms over 5740 records, baseline subtraction
  0.043 ms over 1915.
- Cross-study replication of ranolazine, verapamil and dofetilide: the parallel
  study reads 5.0–9.7 ms higher than the crossover at matched concentration. (Corrected in 0.1.1: dofetilide is a crossover arm and the cause was
  misattributed; see below.)
- ECGRDVQ validated against Vicente 2015 JAHA Table 2 (doi:10.1161/JAHA.114.001615),
  which covers exactly the four bridge drugs: 8 of 8 agreement.
- ECGCIPA validated day-stratified against Vicente 2019 CPT: 6 of 8.
- READING_LIST.md: the papers carrying numbers or data, what each settled, and
  what is still unread.
- 81 tests (`pytest --collect-only -q`; parametrized cases expand, so this is
  larger than the number of `def test_*` functions).
