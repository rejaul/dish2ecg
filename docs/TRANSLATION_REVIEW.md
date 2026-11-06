# From dish to ECG: what differs, what is testable, and what is not

**A literature-grounded review of the sources of difference between an hiPSC-CM field-potential concentration–response and a clinical concentration–QTc relationship, with the limits that public data impose.**

Companion to `docs/RESEARCH_REPORT.md` · 4 November 2026 · literature sweep of 28 PubMed queries (206 records) plus targeted retrieval; every citation resolved against CrossRef. Numbers attributed to this project are reproduced by scripts in this repository; numbers attributed to a paper are from that paper's own abstract or text.

## 1. The question, answered

*Can an in-vitro concentration–response from hiPSC-CM field potentials and a clinical concentration–QTc relationship be passed through the same correction algebra and compared on a common concentration scale — and how far do public data let that be tested?*

**The algebra: yes, and the arithmetic is validated on both sides.** Rate correction, baseline subtraction, control subtraction and a mixed-effects concentration–response model are domain-agnostic operations. Applied to clinical data they reproduce published analyses (per-record Fridericia QTcF matches the ECGCIPA study's own derivation to a median of 0.000 ms, SD 0.190 ms, over 5,740 records; 8 of 8 arm-level agreement with [Vicente 2015]). Applied to in-vitro data the rate-correction step reproduces an independent published implementation (17 of 18 studies within 5 ms of [Millard 2018], median 0.65 ms). The fitted rate exponents on drug-free clinical records, 0.301 and 0.299, sit where Fridericia predicts.

**The concentration scale: yes, and it is the only axis that is genuinely shared.** Multiples of free Cmax is dimensionless and drug-specific, so "1×" means the same thing for every compound on both sides. Critically, the comparison is **invariant to the choice of reference Cmax**, because the same divisor applies to both axes.

**The effect scale: no.** Milliseconds of FPD and milliseconds of QT are different physical quantities (§4.1). The package implements percent-change and standardized-effect conversions, but neither can be used on the public in-vitro corpus because it ships endpoints already baseline- and vehicle-corrected, leaving no baseline level to divide by and no control SD to standardize against. The headline comparison therefore places two incommensurable millisecond scales side by side and says so.

**How far public data let it be tested: four drugs.** That is the entire intersection of the 28-drug in-vitro panel with single-agent clinical arms carrying open concentration–QTc data. At matched free exposure:

| Drug | × free Cmax | Clinic ΔΔQTcF (90% CI) | Dish ΔΔFPDc (90% CI) | Dish ÷ clinic |
|---|---|---|---|---|
| Dofetilide | 1.07 | +73.1 [68.5, 77.7] | +132.0 [93.1, 170.8] | 1.81 |
| Quinidine | 0.29 | +74.9 [69.3, 80.5] | +135.1 [82.2, 187.9] | 1.80 |
| Ranolazine | 0.96 | +11.2 [7.6, 14.8] | +38.8 [22.8, 54.7] | 3.46 |
| Verapamil | 0.50 | +3.5 [−0.1, 7.0] | −86.3 [−98.9, −73.6] | −24.7 |

Three of four agree in sign. The two selective hERG blockers are inflated by almost identical factors (1.81, 1.80), which looks like a property of the assay rather than noise. Verapamil inverts. **Four drugs cannot calibrate anything**; they can only demonstrate that the machinery runs and expose where it breaks.

## 2. Why this is hard: the two systems are not nested

The intuitive model — "the dish is a simplified heart, so it should give a smaller but proportional signal" — is wrong in a specific way. The dish is not a scaled-down heart; it is a *differently wired* preparation measured by a *different physical principle* under a *different exposure*. Differences therefore do not cancel; they compose, and several of them act in opposite directions.

The taxonomy below follows the causal chain from the molecule in the vial to the number in the table. For each source I give the mechanism, what the literature measures, and whether it biases the dish–clinic comparison **up, down, or unknown**.

## 3. Layer 1 — What actually reaches the channel

### 3.1 Nominal is not free, in either system, and not by the same factor

The clinical side corrects for plasma protein binding; the in-vitro side usually does not correct for anything. [Kisitu 2020] states the problem directly: in biological systems drugs are partly bound to biomolecules, only a fraction is free, and because that fraction is set by test-system parameters — protein and lipid content, and **the type of surface material** — comparing nominal concentrations between two different new-approach methods, or between an in-vitro system and human subjects, can lead to faulty conclusions. That is precisely the comparison this package makes.

**Direction: previously assumed, now unknown.** An earlier version of this project argued that free concentration in the well must be below nominal, so correcting it would widen the dish–clinic discrepancy. [Schocken 2018] shows the sign is not established: of 18 drugs that induced arrhythmias or quiescence, **8 did so at lower concentrations in serum-containing medium** than in serum-free — the opposite of what protein binding alone predicts, and the paper calls the trend counterintuitive. Incomplete solubility compounds it: for quinidine and moxifloxacin, measured stock concentrations fell short of nominal.

### 3.2 Serum changes the cells, not just the free fraction

Serum-free medium is used specifically to approximate free drug concentration, but it is not electrophysiologically neutral. [Schocken 2018] measured, in the same cells: interbeat interval 2102 ± 200 ms with serum against 1643 ± 143 ms serum-free, and **baseline FPDc 575 ± 59 ms against 471 ± 44 ms** — a 104 ms baseline difference from medium composition alone, comparable to the largest drug effects being measured. [Wei 2022] reaches the same conclusion from the FDA's own laboratory, attributing the difference partly to abiotic compound loss.

**Direction: unknown, and site-specific.** The ten sites in the in-vitro corpus did not use one medium.

### 3.3 The dish has no liver

A clinical concentration–QTc relationship is driven by whatever circulates, which includes active metabolites. The dish sees only the parent compound. Risperidone's QT effect is predominantly attributable to its 9-hydroxy metabolite paliperidone, and a recent review argues for a metabolite-centric view of QT prolongation generally. The in-vitro assay is structurally blind to this.

**Direction: dish under-reads** for drugs whose metabolites carry the activity.

### 3.4 Which Cmax, and what happens outside the tested range

[Leishman 2020] shows the reference exposure is itself a modelling choice with large consequences: the clinical terfenadine concentration used in CiPA in-silico work was the exposure under metabolic inhibition, **14- to 40-fold above the therapeutic plasma concentration**. This project found the same instability empirically — the CiPA in-silico training table and the Blinova reference differ 13.3-fold for terfenadine, while [Crumb 2016] Table 1 and Blinova agree to a median ratio of 1.000 across 12 shared drugs.

[Chang 2017] bounds where concentration-scaled in-vitro metrics remain trustworthy: uncertainty quantification showed accurate risk stratification by qNet **below 5× Cmax**, but above 10× Cmax drugs could no longer be separated by risk — because for most pharmacology data the measured block was under 60%, preventing reliable IC50 estimation.

**Direction: unknown; magnitude large.** This is why this project evaluates at the *observed* clinical free Cmax rather than a nominal 1×, and flags when a comparison point falls outside the tested in-vitro range (quinidine, 0.288× against a 0.32× floor).

### 3.5 Acute versus chronic exposure

Clinical QT data come from steady state. The in-vitro protocol is typically an acute 30-minute exposure. [Qu 2025] shows this matters for mechanism: FPD prolongation from moxifloxacin appeared within 10 minutes, but detecting the effect of an siRNA targeting hERG required a **48-hour** protocol, because the mechanism is transcript knockdown rather than pore block. Any drug acting through hERG trafficking rather than acute block is invisible on the short protocol.

**Direction: dish under-reads** for trafficking inhibitors and other slow mechanisms.

## 4. Layer 2 — The cell is not an adult ventricular myocyte

### 4.1 Reduced IK1, and what follows from it

This is the single most consequential biological difference and it has a dedicated review. [Goversen 2018] sets it out: unlike adult ventricular cardiomyocytes, iPSC-CMs **beat spontaneously** because of the pacemaker current If and **reduced IK1 density**; in adult cells IK1 finalises repolarisation by stabilising the resting potential while maintaining excitability; the reduced density "contributes to proarrhythmic traits in iPSC-CMs, which leads to an electrophysiological phenotype that might bias drug responses."

Three consequences follow. The resting potential is depolarised, so sodium-channel availability differs. The cells self-pace, so rate is an *outcome* rather than a controlled input (§6.2). And repolarisation leans more heavily on IKr, which is exactly the current most cardiac-safety drugs block.

**Direction: dish over-reads for selective IKr blockers.** The measured 1.81 and 1.80 inflation factors for dofetilide and quinidine are consistent with this, though not proof of it.

### 4.2 Repolarisation reserve is the mechanism, and it has been demonstrated directly

[Braam 2013] is the cleanest demonstration that reserve governs what the assay can see. Eleven reference compounds were applied to hPSC-CM clusters on MEAs; hERG blockers prolonged and calcium blockers shortened FPD as expected, but **two potent IKs blockers had only minor effects despite functional IKs channels** — and became detectable only when the preparation was re-evaluated under conditions mimicking reduced repolarisation reserve. A drug's visibility in the dish depends on the reserve of the preparation, not only on the drug.

[Gaur 2020] develops a quantitative measure of reserve as a proarrhythmia marker, and [Li 2017b] shows the engineering fix: KCNJ2 overexpression (restoring IK1) allowed stable measurement of paced APD prolongation at 0.1 µM E-4031, a concentration that **abolished firing entirely** in untreated cells. The untreated preparation does not merely under-report at that concentration; it stops producing a measurement at all — which is the mechanism behind this project's informative missingness (§6.4).

### 4.3 The preparation is a mixed, drifting population

Three distinct problems, each measured.

**Composition.** Single-cell RNA sequencing of a commercial iCell preparation [Schmid 2020] found it heterogeneous, including a cluster (10.8% of cells) co-expressing cardiac and cell-cycle genes. Monolayers mix ventricular-like, atrial-like and nodal-like cells in proportions that differ by line.

**Line and vendor.** [da Rocha 2020] validated three hiPSC-CM lines in the CiPA Phase II study and reports that each had **distinct baseline electrophysiology — spontaneous beat rate, APD and conduction velocity — and distinct drug responsiveness**. In this project's corpus, cell line is confounded with site and platform across 15 dataset combinations.

**Time.** Maturation continues in culture, so the preparation is a moving target. [Seibertz 2023] tracks electrophysiology and calcium handling through long-term culture. [Guerrelli 2024] quantifies both timescales that matter in practice: within a single recording session, a **22.6% increase in beating rate and 7.7% decrease in FPD over a 20-minute equilibration**; between days 2 and 14 of culture, beating rate −12.7 beats/min, **FPD +257 ms**, spike amplitude +3.3 mV.

**Direction: unknown, large, and the dominant variance component.** In this project the in-vitro 90% intervals are 25–106 ms wide against 7–11 ms clinically, despite 128–278 wells per drug against 21–22 subjects. The variance lives between site × line × platform combinations, not between wells, so adding wells does not help.

## 5. Layer 3 — The tissue that is not there

### 5.1 No transmural gradient, therefore no Tpeak-Tend analogue

The clinical biomarker set splits QT into J-Tpeak (early repolarisation) and Tpeak-Tend (late, dispersion-related). The second has no in-vitro counterpart, because a monolayer has no wall across which to disperse. [Yan 1998] measured the gradient directly in arterially perfused canine wedge preparations: **dispersion of APD90 across the left ventricular wall averaged 51 ± 19 ms at a cycle length of 1000 ms and 64 ± 25 ms at 2000 ms**, with transmural dispersion of *repolarisation time* smaller (34 ± 18 and 45 ± 25 ms) because the endocardial-to-epicardial activation sequence partly offsets it. Tissue resistivity rose 2.8-fold in the M-cell region.

Whether Tpeak-Tend indexes transmural dispersion specifically remains contested in the literature [Antzelevitch 2007], but the structural point is not: a 2D monolayer cannot express a gradient that is defined across a 3D wall, and a field potential has no T-wave morphology to decompose. **The mechanistic discrimination that makes the clinical biomarker set valuable — selective hERG block prolongs J-Tpeak, balanced block does not — has no direct in-vitro expression.**

### 5.2 No activation sequence, no Purkinje system

In the heart, repolarisation time is activation time plus APD, and the His-Purkinje system imposes an activation sequence that partly cancels the APD gradient (the Yan 1998 numbers above). A monolayer activates from wherever it happens to pace. There is no QRS in any meaningful sense, and no analogue of conduction-driven repolarisation heterogeneity.

### 5.3 No autonomic input — and this is not a small term

A dish has no sympathetic or parasympathetic tone. Two human studies quantify what that removes, and the direction is surprising.

[Smith 2007] gave ibutilide to healthy subjects with and without double autonomic blockade. With intact autonomic tone, ibutilide prolonged QTc by **43 ± 8 ms**. Autonomic block alone prolonged QTc by **16.7 ± 9.4 ms**. Ibutilide under double autonomic block produced an **additional 68 ± 5 ms** of prolongation. The paper's conclusion is in its title: autonomic tone *attenuates* drug-induced QT prolongation.

[Diedrich 2002] is consistent: ganglionic blockade with trimethaphan prolonged QTc in normal subjects from 432 ± 6 to 469 ± 7 ms while shortening RR from 995 ± 45 to 670 ± 35 ms.

**Direction: this one makes the dish over-read.** An autonomically denervated preparation is in the state where a hERG blocker's effect is *largest*. This is a mechanistically distinct contributor to the 1.8-fold inflation seen for dofetilide and quinidine, separate from reduced IK1, and I have not seen it invoked in the hiPSC-CM translation literature.

## 6. Layer 4 — The measurement

### 6.1 A field potential is not an action potential and neither is a QT interval

The MEA records an *extracellular* voltage under a syncytium, closer to the first derivative of the transmembrane potential than to the potential itself. [Ernault 2024] addresses exactly this interpretive gap, noting that field potentials resemble extracellular electrograms recorded from myocardium but that **different criteria are conventionally used to interpret the two**, which hampers correct interpretation and translation to the patient. [Hayes 2019] introduced the LEAP method to recover action-potential-like waveforms from intact monolayers, acknowledging the same limitation.

In this project's corpus the problem is compounded: 2 of 15 datasets (15.2% of wells) are voltage-sensing optical rather than MEA, so the column named `ddFPDc` actually holds ΔΔAPD₉₀c for those — two different physical quantities under one name.

### 6.2 Rate and duration are mutually dependent, which breaks the correction's premise

Rate correction assumes rate is an independent variable. In spontaneously beating preparations it is not. [Rast 2016] states the circularity precisely: FPD varies with beating rate, **but beating rate can also be influenced by FPD**, and "this interdependence is not fully understood and therefore mandates careful data analysis and cautious interpretation." A drug that prolongs FPD may thereby slow the rate, which further prolongs the uncorrected FPD.

Pacing removes the circularity. [Wei 2020] and [Patel 2019] both evaluate it; notably, [Patel 2019] found TdP-risk discrimination **identical** for spontaneously beating and 0.8 Hz-paced iPSC-CMs (AUC 0.96, 95% CI [0.9, 1]), and slightly lower at 1 and 1.5 Hz. So pacing fixes a real methodological problem without necessarily improving classification.

**Direction: unknown; a correctness problem rather than a bias.** The public corpus is unpaced, so this project inherits it.

### 6.3 Temperature and protocol

hERG block is temperature-dependent and **not in a consistent direction**: in one survey amiodarone and β-estradiol blocked more potently at 37 °C than at 23 °C (amiodarone IC50 0.30 µM vs 0.56 µM), while ivermectin and frusemide blocked *less* potently at 37 °C (12.52 → 24.41 µM). Ion-channel IC50s measured at room temperature cannot be mapped to a 37 °C assay by a single factor.

### 6.4 Informative missingness: the dish stops reporting exactly where the effect is largest

A well in arrhythmia has no measurable field-potential duration. In the public corpus, of 2,067 missing `ddFPDc` values **1,887 are exactly the wells flagged EAD = 1**, and missingness rises 0.9% → 4.7% → 26.8% → 60.8% across the four concentration indices. This is missingness not at random in the strict sense [Hazewinkel 2022], and it biases complete-case slopes toward zero.

It is also the mechanistic consequence of §4.2: [Li 2017b]'s observation that E-4031 *abolished firing* in IK1-deficient cells is the same phenomenon that produces these blanks.

This project bounds rather than fixes it — every in-vitro fit is reported over all concentrations and over indices 1–3 — and the bound is informative: for the high-risk drugs the two agree closely, but **quinidine's slope changes sign** between the two analyses (−28.3 vs +25.3 ms/decade) because 59.9% of its wells are lost. The *effect* at the comparison point is robust (135.1 vs 119.4 ms); the *slope* is not quotable.

## 7. Layer 5 — The analysis, and the clinical side's own instability

It would be convenient if the clinical anchor were fixed and only the dish were noisy. It is not.

**Between-study.** Ranolazine, verapamil and dofetilide have estimates from two public studies. At matched plasma concentration ECGCIPA reads **5.0 to 9.7 ms higher** than ECGRDVQ every time, and for dofetilide the 90% intervals do not overlap. Decomposition shows the shared placebo group explains only ~2–4 ms; **6–7 ms remains unexplained**.

**Measurement algorithm.** Using PTB-XL+ [Strodthoff 2023], where three algorithms measure the same 21,795 records, two commercial algorithms differ on QT by a mean of −1.87 ms (SD 7.22) but by **more than 10 ms on 21.3% of records**, with a rate-dependent component of −0.0196 ms per ms of RR (≈ −5.9 ms across a 300 ms RR change). The cross-sectional QT/RR exponent fitted on identical records is 0.3355 against 0.3734. A constant bias cancels in ΔΔ; a rate-dependent one does not.

**Correction and model choices.** The two published J-Tpeak corrections differ by up to 21 ms per record (though ≤1.2 ms at the drug level). A random-slope specification versus random-intercept moves individual estimates by up to 4 ms.

**The clinical target is not a fixed point either.** [Wu 2023] reports that **14% of clinical QTc studies for approved and investigational products failed to exclude 10 ms prolongation** — the threshold the whole exercise is calibrated against is itself crossed by a seventh of approved products.

## 8. Layer 6 — Population and endpoint

**Healthy volunteers are the easy case.** The clinical corpus is 104 healthy adults aged 19–35 with no structural heart disease, no electrolyte disturbance and no concomitant QT-prolonging therapy. Metabolic and electrolyte abnormalities are established risk factors in drug-induced long QT [TeBay 2022], and are precisely what the corpus excludes.

**Sex.** Women are at higher risk of drug-induced TdP. [Vicente 2015b] measured sex differences in drug-induced repolarisation changes in the same ECGRDVQ study this project uses, and [Huo 2019] built an hiPSC-CM model system to study the genetic and hormonal components. This project does not stratify by sex; the data support it and it is an unexploited axis.

**The endpoint is a surrogate, in the population least likely to express the outcome.** There are zero arrhythmic events in the clinical corpus and none would be expected. The assay exists to predict torsade; the anchor measures QTc. [Vargas 2021] argues from an industry perspective for a fully integrated nonclinical–clinical assessment rather than the current sequential one, which is a direct response to this mismatch.

**Individual prediction fails, and the dish is not solely at fault.** [Blinova 2019] found no within-person correlation between a subject's own ΔΔAPD₉₀c slope and that subject's clinical ΔΔQTc slope (ρ = −0.106, p = 0.70 dofetilide; ρ = −0.089, p = 0.75 moxifloxacin). The informative control: those subjects' clinical responses to the **two drugs did not correlate with each other** (ρ = 0.096, p = 0.73).

## 9. Summary table — direction of each difference

| # | Source | Direction on dish-vs-clinic | Evidence |
|---|---|---|---|
| 1 | Nominal vs free concentration in well | **Unknown** (was assumed "dish over-potent") | [Kisitu 2020]; [Schocken 2018] 8/18 drugs arrhythmic at *lower* serum concentrations |
| 2 | Serum / medium composition | Unknown, site-specific | [Schocken 2018] 104 ms baseline FPDc shift; [Wei 2022] |
| 3 | No hepatic metabolism | Dish **under**-reads | metabolite-driven QT effects |
| 4 | Reference-Cmax choice | Unknown, up to 14–40× | [Leishman 2020]; [Chang 2017] |
| 5 | Acute vs chronic exposure | Dish **under**-reads slow mechanisms | [Qu 2025] 48 h needed for siRNA |
| 6 | Reduced IK1 / repolarisation reserve | Dish **over**-reads IKr blockers | [Goversen 2018]; [Braam 2013]; [Li 2017b] |
| 7 | Mixed cell population, line, vendor | Unknown; dominant variance | [Schmid 2020]; [da Rocha 2020] |
| 8 | Time in culture / within-session drift | Unknown; large | [Guerrelli 2024] FPD +257 ms over 12 days |
| 9 | No transmural gradient | No Tpeak-Tend analogue at all | [Yan 1998] 51 ± 19 ms APD90 dispersion |
| 10 | No activation sequence / Purkinje | No QRS analogue | [Yan 1998] |
| 11 | **No autonomic tone** | Dish **over**-reads | [Smith 2007] tone attenuates by ~25 ms; [Diedrich 2002] |
| 12 | Field potential ≠ AP ≠ QT | Not commensurable | [Ernault 2024]; [Hayes 2019] |
| 13 | Rate–duration circularity | Correctness problem | [Rast 2016]; [Patel 2019] |
| 14 | Temperature | Drug-specific, both directions | IC50 shifts at 23 vs 37 °C |
| 15 | Informative missingness | Complete-case slopes → zero | this project: 1,887/2,067 nulls are EAD wells |
| 16 | Clinical between-study instability | ±5–10 ms, mostly unexplained | this project |
| 17 | Interval-extraction algorithm | Rate-dependent, ~6 ms | this project, PTB-XL+ |
| 18 | Healthy-volunteer population | Dish-clinic pair both optimistic | [TeBay 2022] |
| 19 | Sex | Unmodelled | [Vicente 2015b]; [Huo 2019] |
| 20 | QTc surrogate, not TdP | Structural | [Vargas 2021]; [Wu 2023] 14% of studies fail 10 ms |

**The honest reading of this table: at least two mechanisms (6 and 11) push the dish to over-read selective IKr blockers, which is the direction actually observed (1.81×, 1.80×). But at least three (1, 2, 7) have unknown sign and plausibly comparable magnitude. The 1.8-fold factor should not be treated as a calibration constant.**

## 10. Open problems, ranked by what would move the field

1. **No public well-level dataset carries a measured free concentration.** Every in-vitro concentration axis in the field is nominal. Until someone publishes measured free concentrations per well per site, §3.1 cannot be closed and the shared concentration axis is shared only by assumption. *This is the single highest-value missing measurement.*
2. **The four-drug ceiling.** Widening the bridge needs single-agent clinical arms with open concentration–QTc data for drugs in the in-vitro panel. None exist beyond the four.
3. **The 6–7 ms unexplained between-study shift**, which bounds how precisely any dish-to-clinic relationship can be known from one clinical study.
4. **Autonomic tone as a translation factor** (§5.3). The human measurements exist [Smith 2007]; nobody appears to have used them to *correct* an in-vitro over-read. A testable hypothesis: part of the 1.8× inflation is denervation, not immaturity.
5. **Informative missingness has no principled estimator here.** A tipping-point analysis — assign missing wells the observed maximum plus δ and report the δ at which conclusions change — is specified but not implemented.
6. **Tpeak-Tend has no in-vitro analogue**, so the clinically most discriminating biomarker cannot cross the bridge. Engineered 3D tissue with a transmural gradient is the obvious route and is outside what public data can test.
7. **Reserve-matched comparison.** [Braam 2013] implies the right comparison is not dish-vs-clinic at matched concentration but at matched repolarisation reserve. No public dataset reports reserve per well.
8. **Sex stratification** is supported by the public clinical data and not done.
9. **The endpoint gap.** Zero TdP events in the corpus; QTc is all there is to anchor on.

## 11. What this review does not establish

- The literature sweep was 28 PubMed queries plus targeted retrieval, in November 2026. It is not a systematic review: no PRISMA protocol, no dual screening, no formal inclusion criteria, and no assessment of study quality. Absence of a topic here is weak evidence of absence in the literature.
- Claims attributed to papers come from abstracts and, where retrieved, full text. Where only an abstract was read, the numbers quoted are those the authors chose to put in the abstract.
- The direction assignments in §9 are my reading of the mechanisms, not measured decompositions. Only sources 15–17 are quantified within this project.
- No new experiment is proposed that this project could run on public data; items 1, 4, 6 and 7 all require data that does not exist publicly.

## References

- **[Antzelevitch 2007]** Antzelevitch C et al. (2007) Does Tpeak-Tend provide an index of transmural dispersion of repolarization?. *Heart Rhythm* 4:1114-1116. doi:10.1016/j.hrthm.2007.05.028
- **[Blinova 2019]** Blinova K et al. (2019) Clinical Trial in a Dish: Personalized Stem Cell-Derived Cardiomyocyte Assay Compared With Clinical Trial Results for Two QT-Prolonging Drugs. *Clinical and Translational Science* 12:687-697. doi:10.1111/cts.12674
- **[Braam 2013]** Braam SR et al. (2013) Repolarization reserve determines drug responses in human pluripotent stem cell derived cardiomyocytes. *Stem Cell Research* 10:48-56. doi:10.1016/j.scr.2012.08.007
- **[Chang 2017]** Chang KC et al. (2017) Uncertainty Quantification Reveals the Importance of Data Variability and Experimental Design Considerations for in Silico Proarrhythmia Risk Assessment. *Frontiers in Physiology* 8. doi:10.3389/fphys.2017.00917
- **[Crumb 2016]** Crumb WJ et al. (2016) An evaluation of 30 clinical drugs against the comprehensive in vitro proarrhythmia assay (CiPA) proposed ion channel panel. *Journal of Pharmacological and Toxicological Methods* 81:251-262. doi:10.1016/j.vascn.2016.03.009
- **[Diedrich 2002]** Diedrich A et al. (2002) Modulation of QT Interval During Autonomic Nervous System Blockade in Humans. *Circulation* 106:2238-2243. doi:10.1161/01.cir.0000035241.76918.6c
- **[Ernault 2024]** Ernault AC et al. (2024) Interpretation of field and LEAP potentials recorded from cardiomyocyte monolayers. *American Journal of Physiology-Heart and Circulatory Physiology* 326:H800-H811. doi:10.1152/ajpheart.00463.2023
- **[Gaur 2020]** Gaur N et al. (2020) Validation of quantitative measure of repolarization reserve as a novel marker of drug induced proarrhythmia. *Journal of Molecular and Cellular Cardiology* 145:122-132. doi:10.1016/j.yjmcc.2020.04.019
- **[Goversen 2018]** Goversen B et al. (2018) The immature electrophysiological phenotype of iPSC-CMs still hampers in vitro drug screening: Special focus on I K1. *Pharmacology & Therapeutics* 183:127-136. doi:10.1016/j.pharmthera.2017.10.001
- **[Guerrelli 2024]** Guerrelli D et al. (2024) hiPSC-CM electrophysiology: impact of temporal changes and study parameters on experimental reproducibility. *American Journal of Physiology-Heart and Circulatory Physiology* 327:H12-H27. doi:10.1152/ajpheart.00631.2023
- **[Hayes 2019]** Hayes HB et al. (2019) Novel method for action potential measurements from intact cardiac monolayers with multiwell microelectrode array technology. *Scientific Reports* 9. doi:10.1038/s41598-019-48174-5
- **[Hazewinkel 2022]** Hazewinkel A et al. (2022) Sensitivity to missing not at random dropout in clinical trials: Use and interpretation of the trimmed means estimator. *Statistics in Medicine* 41:1462-1481. doi:10.1002/sim.9299
- **[Huo 2019]** Huo J et al. (2018) Sex-Related Differences in Drug-Induced QT Prolongation and Torsades de Pointes: A New Model System with Human iPSC-CMs. *Toxicological Sciences*. doi:10.1093/toxsci/kfy239
- **[Kisitu 2020]** Kisitu J (2020) Chemical concentrations in cell culture compartments (C5) - free concentrations. *ALTEX*. doi:10.14573/altex.2008251
- **[Leishman 2020]** Leishman DJ (2020) Improving prediction of torsadogenic risk in the CiPA in silico model by appropriately accounting for clinical exposure. *Journal of Pharmacological and Toxicological Methods* 101:106654. doi:10.1016/j.vascn.2019.106654
- **[Li 2017b]** Li M et al. (2017) Overexpression of KCNJ2 in induced pluripotent stem cell-derived cardiomyocytes for the assessment of QT-prolonging drugs. *Journal of Pharmacological Sciences* 134:75-85. doi:10.1016/j.jphs.2017.05.004
- **[Millard 2018]** Millard D et al. (2018) Cross-Site Reliability of Human Induced Pluripotent stem cell-derived Cardiomyocyte Based Safety Assays Using Microelectrode Arrays: Results from a Blinded CiPA Pilot Study. *Toxicological Sciences* 164:550-562. doi:10.1093/toxsci/kfy110
- **[Patel 2019]** Patel D et al. (2019) Assessment of Proarrhythmic Potential of Drugs in Optogenetically Paced Induced Pluripotent Stem Cell-Derived Cardiomyocytes. *Toxicological Sciences* 170:167-179. doi:10.1093/toxsci/kfz076
- **[Qu 2025]** Qu Y et al. (2025) Delayed Repolarization Caused by                     hERG                     Block With Different Drug Modalities Can Be Detected in Stem Cell-Derived Cardiomyocytes: Incubation Time Matters. *Clinical and Translational Science* 18. doi:10.1111/cts.70283
- **[Rast 2016]** Rast G et al. (2016) Influence of field potential duration on spontaneous beating rate of human induced pluripotent stem cell-derived cardiomyocytes: Implications for data analysis and test system selection. *Journal of Pharmacological and Toxicological Methods* 82:74-82. doi:10.1016/j.vascn.2016.08.002
- **[Schmid 2020]** Schmid C et al. (2020) Characterization of iCell cardiomyocytes using single-cell RNA-sequencing methods. *Journal of Pharmacological and Toxicological Methods* 106:106915. doi:10.1016/j.vascn.2020.106915
- **[Schocken 2018]** Schocken D et al. (2018) Comparative analysis of media effects on human induced pluripotent stem cell-derived cardiomyocytes in proarrhythmia risk assessment. *Journal of Pharmacological and Toxicological Methods* 90:39-47. doi:10.1016/j.vascn.2017.11.002
- **[Seibertz 2023]** Seibertz F et al. (2023) Electrophysiological and calcium-handling development during long-term culture of human-induced pluripotent stem cell-derived cardiomyocytes. *Basic Research in Cardiology* 118. doi:10.1007/s00395-022-00973-0
- **[Smith 2007]** SMITH AH et al. (2007) Autonomic Tone Attenuates Drug-Induced QT Prolongation. *Journal of Cardiovascular Electrophysiology* 18:960-964. doi:10.1111/j.1540-8167.2007.00901.x
- **[Strodthoff 2023]** Strodthoff N et al. (2023) PTB-XL+, a comprehensive electrocardiographic feature dataset. *Scientific Data* 10. doi:10.1038/s41597-023-02153-8
- **[TeBay 2022]** TeBay C et al. (2022) Metabolic and electrolyte abnormalities as risk factors in drug-induced long QT syndrome. *Biophysical Reviews* 14:353-367. doi:10.1007/s12551-022-00929-7
- **[Vargas 2021]** Vargas HM et al. (2020) Time for a Fully Integrated Nonclinical-Clinical Risk Assessment to Streamline QT Prolongation Liability Determinations: A Pharma Industry Perspective. *Clinical Pharmacology & Therapeutics* 109:310-318. doi:10.1002/cpt.2029
- **[Vicente 2015]** Vicente J et al. (2015) Comprehensive T wave Morphology Assessment in a Randomized Clinical Study of Dofetilide, Quinidine, Ranolazine, and Verapamil. *Journal of the American Heart Association* 4. doi:10.1161/JAHA.114.001615
- **[Vicente 2015b]** Vicente J et al. (2015) Sex differences in drug-induced changes in ventricular repolarization. *Journal of Electrocardiology* 48:1081-1087. doi:10.1016/j.jelectrocard.2015.08.004
- **[Wei 2020]** Wei F et al. (2020) Effects of Electrical Stimulation on hiPSC-CM Responses to Classic Ion Channel Blockers. *Toxicological Sciences* 174:254-265. doi:10.1093/toxsci/kfaa010
- **[Wei 2022]** Wei F et al. (2022) Effects of Serum and Compound Preparation Methods on Delayed Repolarization Evaluation With Human iPSC-CMs. *Toxicological Sciences* 188:48-61. doi:10.1093/toxsci/kfac043
- **[Wu 2023]** Wu WW et al. (2023) ICH S7B In Vitro Assays Do Not Address Mechanisms of QTC Prolongation for Peptides and Proteins - Data in Support of Not Needing Dedicated QTC Studies. *Clinical Pharmacology & Therapeutics* 114:1332-1341. doi:10.1002/cpt.3047
- **[Yan 1998]** Yan GX et al. (1998) Characteristics and Distribution of M Cells in Arterially Perfused Canine Left Ventricular Wedge Preparations. *Circulation* 98:1921-1927. doi:10.1161/01.cir.98.18.1921
- **[da Rocha 2020]** da Rocha AM et al. (2019) Detection of Drug-Induced Torsades de Pointes Arrhythmia Mechanisms Using hiPSC-CM Syncytial Monolayers in a High-Throughput Screening Voltage Sensitive Dye Assay. *Toxicological Sciences* 173:402-415. doi:10.1093/toxsci/kfz235

*Reference list contains only works cited above; the full project bibliography is in REFERENCES.md.*
