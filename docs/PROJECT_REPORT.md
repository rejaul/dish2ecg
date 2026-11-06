# From the dish to the ECG

*A public-data comparison of stem-cell heart cells and clinical QT results*

Project report · dish2ecg version 0.2.2 · 5 November 2026

## Summary

*Readers who want the physiology before the methods should start at section 12, and those who want the reading list should start at section 17.*

Safety teams now screen new drugs for torsades risk partly in a dish, using heart cells grown from human induced pluripotent stem cells (hiPSC-CMs). The dish gives a field potential duration (FPD). It looks like the QT interval on a patient's ECG, but it is not one, and published work often puts the two side by side without saying how each was corrected or on what concentration scale. We wanted to know how far that comparison can be tested with public data alone.

We wrote an open Python package, dish2ecg, that sends dish data and clinical ECG data through the same correction steps and compares them as multiples of the free drug concentration that patients reach. Before relying on each stage we checked it against a published analysis. The clinical side matched published tables in 8 of 8, 5 of 6 and 6 of 8 comparisons, and our dish-side heart-rate correction fell within 5 ms of an independent implementation in 17 of 18 studies.

Only four drugs have both dish data and clinical concentration data: dofetilide, quinidine, ranolazine and verapamil. At the same free concentration the dish gives the right direction for three of them and overstates the effect by a factor of 1.8 to 3.5. For verapamil it goes the wrong way: +3.5 ms in volunteers and −86 ms in the dish. Sorting 28 drugs by torsades risk is a more forgiving task. On the 16 drugs that the FDA's in-silico score (qNet) was not fitted on, the dish ranks drugs about as well as qNet does (Spearman correlation with risk class 0.81 against 0.84) and is as good at clearing low-risk drugs. qNet is better at flagging the high-risk ones.

Four problems stay open. Public data cannot separate the type of dish measurement from the laboratory and the culture medium. The clinical recordings are only 10 seconds long, so the usual correction for QT "memory" cannot be applied. Two clinical studies disagree by 6 to 7 ms for the same drugs and we cannot explain it. And a four-drug bridge cannot support a calibration. Along the way we withdrew four of our own claims after testing them; they are listed in section 5.

## 1. The question

A cardiac safety package asks one thing twice. Does this compound delay repolarisation of the ventricle, and by how much, at the concentration a patient will actually have in the blood? In the clinic the answer comes from the QT interval of a 12-lead ECG, corrected for heart rate (QTc). In the laboratory, since the CiPA initiative (Sager 2014; Colatsky 2016; Fermini 2016), part of the answer comes from a plate of hiPSC-CMs sitting on electrodes. That plate records a field potential whose duration (FPD) is corrected for beat rate in the same spirit.

The two readouts are routinely placed side by side. A reader can rarely check how each was corrected, which concentration scale was used, or whether the dish number is being read as if it were milliseconds of QT. The largest public dish dataset is the CiPA myocyte study of Blinova and colleagues (2018): 28 drugs tested at 10 sites on two commercial cell lines. The largest public source of clinical concentration–QT data is PhysioNet, which holds three early-phase studies from the FDA and its collaborators.

Our question was whether both can be passed through one correction algebra and compared on one concentration scale, and how far public data let us test that. We did not set out to show that the dish predicts QT. At the start we made one rule that the code enforces: a dish endpoint is never expressed as a QT interval in milliseconds, and the package raises an error if asked to do so.

*Figure 1. What the two dish readouts measure. (a) The field potential from an electrode array: FPD runs from the depolarisation spike to the peak of the repolarisation wave, which mainly reflects the potassium current IKr. (b) The action potential from a voltage-sensitive optical probe: APD90 is the time to 90 % repolarisation. Neither is a QT interval.* (figure in the HTML edition of this report)

## 2. Data

Everything came from public sources. We did not redistribute any third-party data; the repository holds code, derived results and pointers to the originals.

| Source | What it holds | Size |
|---|---|---|
| CiPA myocyte study, Blinova et al. 2018 (workbook supplied by the project owner) | Well-level corrected change in field potential duration (ddFPDc), arrhythmia flag, 4 concentrations per drug | 28 drugs, 8,859 wells, 10 sites, 6 platforms, 2 cell lines (iCell², Cor.4U), 15 dataset-by-platform combinations |
| ECGRDVQ (PhysioNet) | Dofetilide, quinidine, ranolazine, verapamil and placebo: intervals, plasma levels, waveforms | 22 subjects, crossover |
| ECGDMMLD (PhysioNet) | Dofetilide with or without mexiletine or lidocaine; moxifloxacin with or without diltiazem | 22 subjects, crossover |
| ECGCIPA (PhysioNet) | Ranolazine, verapamil, lopinavir/ritonavir, chloroquine (parallel arms); dofetilide with or without diltiazem (crossover) | 60 subjects |
| Millard et al. 2018, Supplementary Table 1 | Beat period, FPD and corrected FPD for 18 studies | 18 rows |
| FDA CiPA repository (pinned commit 539ff80) | Ion-channel potency fits for 28 drugs; 2,000 bootstrap hERG fits per drug; the ORd-CiPA model | 28 drugs |
| Crumb et al. 2016 Table 1, product labels | Free Cmax; plasma protein binding | 30 drugs; 15 reference rows |

The clinical anchor is 11 molecules in 104 healthy volunteers. Nobody in these studies had an arrhythmia, which matters later. In the dish table, 2,067 of 8,859 values are missing, and 1,887 of those are exactly the wells flagged as early afterdepolarisations (EADs), that is, wells that stopped beating normally.

## 3. Methods

### 3.1 One cascade for both sides

The core idea is plain on purpose. A raw interval is corrected for rate, then the subject's or well's own baseline is subtracted (Δ), then the placebo or vehicle change is subtracted (ΔΔ), then a concentration–response model is fitted, and the result is put on a shared concentration scale. Each step is recorded in a provenance tuple attached to the number. The code refuses steps in the wrong order, refuses to subtract a control from something that has no baseline, and refuses to subtract twice.

![Figure 2. The processing pipeline, drawn. Rows are the shared cascade steps; columns are the three arms. Every panel below the top row is this repository's own data: the rate-correction row shows each of the 18 Millard studies moving from FPD to FPDc and 439 drug-free clinical records; the subtraction row shows the real dofetilide waterfall at 2.5 h; the model row shows the three fitted concentration-responses with 90 % bands, with the ICH E14 10 ms line drawn on the clinical panel and the in-silico curve cut off where the mid-myocardial cell stops repolarising. The top row is schematic, since no raw waveform is redistributed here. Red marks the places the method refuses to go: expressing FPD as a QT interval, and treating a nominal dish concentration as a measured free one.](../results/pipeline_diagram.png)

### 3.2 Rate correction

We used the power form Ic = I / RR^α, with RR in seconds. Fridericia's α = 1/3 is the default. Exponents fitted on drug-free records (baseline and placebo only, so that the drug effect is not absorbed) came out at 0.301 in ECGRDVQ and 0.299 in ECGDMMLD, which is what healthy volunteers should give and acts as a check on the fitting code. For J-Tpeak we used the published exponent 0.58 (Johannesen 2014) and confirmed it on 5,740 ECGCIPA records, where a regression recovers 0.5799 with a residual SD of 0.15 ms. This exponent is not ours.

### 3.3 Concentration–response models

On the clinical side we fitted a linear mixed model of ΔΔQTcF against concentration with a random intercept per subject, and read the effect at the geometric-mean observed Cmax with a two-sided 90 % interval. Over 200 simulated crossover studies the interval covered the true effect 86.5 % of the time against a nominal 90 %, so our clinical intervals are slightly too narrow.

On the dish side we used the same model with the dataset (site, cell line and platform) as the random effect and log10 of the multiple of free Cmax as the predictor. A log scale was needed because tested concentrations span up to 3,500-fold and a linear fit is dominated by its highest point. Wells that stop beating cannot be measured, and the loss rises steeply with dose (0.9 %, 4.7 %, 26.8 % and 60.8 % across the four concentrations). The fits are therefore reported twice, once over all four concentrations and once over the lowest three. We bound the problem; we do not claim to fix it.

### 3.4 The bridge between dish and clinic

Plasma concentration in ng/mL was converted to nM with the molecular weight, then multiplied by the fraction unbound to get a free concentration, then divided by the free Cmax to get multiples of Cmax. Both sides were evaluated at the geometric-mean free Cmax that the volunteers actually reached, not at a nominal 1×. Quinidine only reached 0.29× and verapamil 0.50× of the reference value, so a nominal 1× would have meant extrapolating the clinical fit. Because the same reference scales both axes, the bridge does not depend on which published Cmax table is used. One asymmetry survives all of this and it is not small: the clinical axis is a measured free plasma concentration, while the dish axis is the nominal concentration added to the well. Section 7 takes that up.

### 3.5 Waveforms

The PhysioNet studies turned out to include raw waveforms and median beats, not only interval tables. We decoded them (WFDB formats 16 and 32, per-microvolt gains), detected R-peaks, and re-measured the T-wave end by a tangent rule. Within each 10-second strip we computed the spread of QTcF between the slowest and the fastest beat.

### 3.6 In-silico work

We ran the O'Hara–Rudy human ventricular model (O'Hara 2011) in the myokit package: 408 simulations covering endocardial, epicardial and mid-myocardial cells, each of the four bridge drugs and 34 concentrations per drug, with no failures. For the risk comparison we used the FDA's ORd-CiPA v1.0 model (Dutta 2017) and computed qNet for 28 drugs, each at four doses and for 20 bootstrap hERG parameter sets: 2,240 simulations, no failures. A coupled one-dimensional cable could not be run, for a software reason described in section 7.

### 3.7 Checking ourselves

The repository has 133 automatic tests (counted with pytest's collect-only mode, not from memory). Every result file has a producing script or is labelled as made by hand. An audit near the end found 14 of 38 result files with no producer; nine are now reproducible, one was deleted and four are labelled ad hoc in `results/README.md`.

## 4. Results

### 4.1 The code reproduces published analyses

| Stage | Comparison | Result |
|---|---|---|
| Source data | Blinova workbook against all 420 published panels | r = 0.999993; median absolute difference 0.30 ms |
| Clinical cascade | Vicente 2015, Table 2 | 8 of 8 |
| Clinical cascade | Vicente 2019, day-stratified | 6 of 8 (both misses are verapamil, in opposite directions) |
| Clinical cascade | Vicente 2016, Table 3 | 5 of 6 |
| Clinical cascade | Johannesen 2015, dofetilide alone | 39.0 and 22.2 ms against their 39.1 and 22.4 ms |
| Per-record algebra | QTcF recomputed on 5,740 ECGCIPA records | median difference 0.000 ms, SD 0.19 ms |
| Dish rate correction | Millard 2018, 18 studies | 17 of 18 within 5 ms; median 0.65 ms |

Two points deserve a sentence. The dofetilide line closed what had been our one unexplained gap: it was a disagreement between two published analyses of the same trial, and our value matches the primary paper. In the dish rate-correction check, study 9 is off by 22.0 ms and we do not know why; it stays visible as a test that asserts the discrepancy.

### 4.2 The four-drug bridge

![Figure 3. Change in repolarisation duration against multiples of free Cmax for the four drugs that exist in both data sets. Clinical ΔΔQTcF in blue, dish ΔΔFPDc in orange, matched-exposure readout in each panel. The dish agrees in sign for three drugs and reverses for verapamil.](../results/dish_vs_clinic.png)

| Drug | Free Cmax multiple | Clinic ΔΔQTcF (ms) | Dish ΔΔFPDc (ms) | Dish ÷ clinic |
|---|---|---|---|---|
| Dofetilide | 1.07 | +73.1 [68.5, 77.7] | +132.0 [93.1, 170.8] | 1.81 |
| Quinidine | 0.29 | +74.9 [69.3, 80.5] | +135.1 [82.2, 187.9] | 1.80 |
| Ranolazine | 0.96 | +11.2 [7.6, 14.8] | +38.8 [22.8, 54.7] | 3.46 |
| Verapamil | 0.50 | +3.5 [−0.1, 7.0] | −86.3 [−98.9, −73.6] | −24.7 |

The two pure hERG blockers are inflated by almost the same factor (1.81 and 1.80), which looks like a property of the assay and not of the drugs. The dish intervals are several times wider than the clinical ones despite 128 to 278 wells per drug. The reason is the variance between sites, platforms and cell lines, not noise between wells. Four drugs is the whole overlap; a search of PhysioNet, Zenodo and the Scientific Data descriptors found no other public well-level field-potential dataset, and we found no other public single-agent concentration–QTc set.

### 4.3 Why verapamil reverses

At the Crumb free Cmax, dofetilide blocks hERG by 24 % and has little effect on the calcium and sodium currents. Verapamil blocks hERG by 15 % and the L-type calcium current by 16 %. The usual reading is that in volunteers the calcium block offsets the potassium block, so QTc barely moves. A likely reason for the dish result is that immature hiPSC-CMs have little of the inward rectifier IK1 and beat spontaneously, so the calcium current dominates their plateau and blocking it shortens the field potential. We did not test this.

We tested whether verapamil's own channel profile is enough to cause the reversal by running it through the ORd model. It is not: the model prolongs the action potential slightly for verapamil, as the clinic does, and gets the sign right for all four drugs. So the dish is the part that is off. The model is far from perfect on size. Ranolazine comes out at about +37 ms in the mid-myocardial cell against +11 ms in volunteers, close to the dish's +39 ms. We read the in-silico arm as a check on direction and not as a prediction of QTc.

The uncoupled model also says something about the QT surrogate. The mid-myocardial (M-cell) action potential duration was the longest in 121 of 136 conditions; the 15 exceptions all sit at 3.78× free Cmax or above. So an uncoupled "QT" is simply the longest cell, with an implied dispersion of about 102 ms where a real wall has 30 to 40 ms. The common mapping that takes a weighted average of action potential durations has no mechanism behind it. The control durations were 265.5 ms (endocardial), 228.0 ms (epicardial) and 329.9 ms (mid).

### 4.4 Endpoint, laboratory and medium cannot be separated

Blinova 2018 stores two different measurements in one column, `ddFPDc`. Thirteen of the 15 datasets are electrode FPD. The other two are the same optical platform on the two cell lines (15.2 % of wells): APD90 measured with a voltage-sensitive probe. The paper calls the two endpoints "equivalent measures". They are not equivalent as a scale. The optical arm reads 1.44 times the consensus of the other nine sites (interval 1.15 to 1.72), the highest of the ten. In the paired comparison, optical = 1.507 × electrode + 24.8 ms (n = 97).

The difficulty is that the optical platform was used at one site only, and that site is also the only one that used serum-free medium. Schocken and colleagues (2018) showed that medium alone moved baseline corrected FPD by 104 ms in iPSC-CMs (575 against 471 ms). So endpoint, site and medium move together, and no public dataset separates them. We can say that the optical arm reads higher; we cannot say why. This is a statement about what the data can and cannot identify, and we think it is more useful than a number.

*Figure 4. Three rungs from cell to patient. Left: per-site gain relative to the other nine sites, with the single optical site highest. Middle: optical against electrode readings of the same drug and concentration, slope 1.51. Right: clinical, dish and in-silico effects for the four bridge drugs.* (figure in the HTML edition of this report)

### 4.5 The clinical side has its own floor

Each clinical ECG record in the waveform data is 10 seconds long. The QT interval takes one to two minutes to settle after a change in heart rate, and 10 seconds is only 8 to 17 % of one adaptation time constant. Most of the RR history that sets the QT therefore lies before the strip begins. This makes the hysteresis correction unidentifiable here, not merely unimplemented.

The cost is measurable. Across 238 records, QTcF computed on the slowest and the fastest beat within one strip differs by a median of 12.6 ms, and by more than 10 ms in 66.8 % of records. The studies take replicate ECGs at each time point, which reduces the effect on a mean. Even so, it shows how little a single 10-second strip fixes about steady-state QT.

The T-wave end is also rate dependent. Our tangent rule against the study's own markings differs by a median of −4.5 ms in ECGRDVQ and −7.1 ms in ECGCIPA, and the difference changes with RR (−0.032 ms per ms of RR, p = 0.026). A constant offset between two delineation methods cancels in a double-delta analysis. A rate-dependent offset does not, if heart rate differs between baseline and treatment. PTB-XL+, where two commercial algorithms measured the same 21,795 ECGs, shows the same sign: a difference of about −5.9 ms across a 300 ms change in RR.

![Figure 5. Waveform findings from 238 clinical records. Left: beat-to-beat RR inside single ten-second recordings. Middle: spread of QTcF between each window's slowest and fastest beat. Right: difference between our tangent T-end and the study's, against RR.](../results/waveform_findings.png)

### 4.6 Sorting 28 drugs by risk

The assay exists to sort drugs into torsades risk classes (8 high, 11 intermediate, 9 low in the CiPA coding). We compared dish endpoints against the in-silico qNet score on all 28 drugs, of which 12 were used to train qNet and 16 were held out.

We first tried this on the 12 training drugs only and concluded that the dish "detects danger but does not grade safety". That conclusion was wrong. It came from looking at the drugs on which qNet had been tuned, with only 12 points. Once the FDA's bootstrap inputs for all 28 drugs turned out to be public, we could test the 16 held-out drugs. Our qNet reproduces the published point fit on 11 of the 12 training drugs to within 0.0006 µC/µF.

| On the 16 held-out drugs | Spearman ρ with risk class | AUC, high vs rest | AUC, low vs rest |
|---|---|---|---|
| qNet (in silico) | 0.840 [0.52, 0.95] | 1.000 | 0.909 |
| ddFPDc at Cmax (dish) | 0.809 [0.55, 0.93] | 0.917 | 0.945 |
| max ddFPDc (dish) | 0.722 [0.38, 0.90] | 0.812 | 0.945 |
| Arrhythmia-event fraction (dish) | 0.717 [0.37, 0.90] | 0.875 | 0.891 |

qNet separates every high-risk drug from every other drug, on all 28. That is the only comparison in the table with no overlap. On ranking and on clearing low-risk drugs the dish and qNet cannot be told apart at n = 16: the intervals overlap almost completely. qNet's two false alarms, both held out, are metoprolol (0.0702) and tamoxifen (0.0781).

*Figure 6. (a) qNet for 28 drugs, median over 20 bootstrap hERG fits; open circles are the 16 held-out drugs. Every high-risk drug falls on the left of the dashed line. (b) Areas under the ROC curve on the 16 held-out drugs.* (figure in the HTML edition of this report)

One caution on reading the table. Two of qNet's traps cost us time and may cost someone else. The current set that qNet integrates excludes the sodium current, and the L-type calcium current acts as the calcium carrier only; including INa gave 0.0044 where the correct value is 0.0785. The trapping constant Kt = 3.5e-5 is fixed in CiPA but defaults to zero in the model file, which switches trapping off and made dofetilide at 4× prolong by 28 ms instead of 130 ms.

### 4.7 A 6 to 7 ms shift we cannot explain

Ranolazine, verapamil and dofetilide were measured in both ECGRDVQ and ECGCIPA. At a matched plasma concentration ECGCIPA reads higher every time: by 8.2 ms for ranolazine, 5.0 for verapamil and 9.7 for dofetilide. ECGCIPA subtracts one shared placebo group from every arm, which accounts for roughly 2 to 4 ms and for essentially all of the verapamil shift. For ranolazine and dofetilide 6 to 7 ms remains; with no placebo subtraction at all the shifts are 6.0 and 7.0 ms.

We tested the most obvious explanation, the T-end rule, on 3,163 median beats from ECGRDVQ and 2,515 from ECGCIPA. It moves the drug effect by −0.94 and −1.83 ms in ECGRDVQ (ranolazine and dofetilide) and by −0.61 and +0.25 ms in ECGCIPA. The difference between the studies is −0.32 ms for ranolazine and −2.08 ms for dofetilide, at most about a third of the dofetilide shift and a twentieth of the ranolazine shift. The placebo-corrected change in RR is only −33 to +18 ms, so a rate-dependent rule would predict a movement of just +0.3 to +1.0 ms. So the shift stays open. We had briefly written it up as a live contributor; that claim is withdrawn.

### 4.8 Free Cmax and fraction unbound for the reference drugs

The bridge needs the fraction unbound (fu) for each drug. For the four bridge drugs it is hard-coded from the literature: dofetilide 0.35, quinidine 0.16, ranolazine 0.38, verapamil 0.10. For four other CiPA drugs that had no fu in our table we looked it up late in the project:

| Drug | fu (point) | Range | Basis |
|---|---|---|---|
| Bepridil | 0.005 | 0.003–0.010 | Label: more than 99 % bound; equilibrium dialysis (Pritchard 1985): more than 99.7 % bound; binding saturates |
| Chlorpromazine | 0.055 | 0.03–0.08 | 92–97 % bound (drug monograph, not an FDA label) |
| Cisapride | 0.0225 | 0.020–0.025 | Label: 97.5–98 % bound |
| Terfenadine | 0.03 | none given | Label (1995): 97 % bound |

None of the four has a current FDA label, so the text comes from reproductions on drug-information sites. Bepridil's point value is a bounded estimate and not a measurement, because its binding depends on concentration. These four drugs have no patient data and are not in the bridge, so no result above depends on them.

The search turned up one finding worth keeping. Terfenadine's free Cmax is 4 nM in the CiPA table and 0.3 nM in Crumb and Blinova, a 13-fold difference. With fu = 0.03 these imply total concentrations of about 63 and 4.7 ng/mL, against a reported normal-dose peak near 1.5 ng/mL. A Redfern 2003 excerpt, which we read only as a search snippet and not in full, says the 4 nM value is a concentration reached with a P450 inhibitor present. If so the two tables describe different exposures, and a "multiples of Cmax" axis means something different for terfenadine depending on which one is used. We also corrected a statement in an earlier note: Crumb 2016 and Li 2017 do not tabulate protein binding.

## 5. What we got wrong

Four claims were made and then withdrawn after a test went against them. We list them because the record of corrections is part of the result.

1. We said that correcting for free drug concentration would widen the gap between dish and clinic. Schocken 2018 shows serum-containing medium made 8 of 18 drugs arrhythmic at lower concentrations, the opposite of what protein binding alone predicts. The sign is not established.
2. We blamed unpaired placebo subtraction in a parallel design for the difference between two clinical studies. That did not survive a direct test.
3. We said the dish detects danger but does not grade safety. This was an artefact of using only the 12 training drugs.
4. We said the T-end rule was a live contributor to the 6 to 7 ms shift. The beat-level test showed it is not.

Smaller errors were also fixed: a rate-correction median quoted as 0.6 ms that is 0.65 ms; a paper credited with naming one "best" predictor when it names three; a test count that no longer matched the suite; and an earlier note claiming that two papers tabulate plasma binding when they do not.

## 6. Discussion

The dish and the clinic agree on direction for most drugs and disagree on size. That is not a surprise. A denervated, immature, single-layer preparation with no His–Purkinje system should over-read a hERG blocker, and a human study found that autonomic tone attenuates drug-induced QT prolongation: ibutilide prolonged QTc by 43 ± 8 ms with intact tone and by 68 ± 5 ms, on top of the blockade's own effect, under double autonomic block (Smith 2007). That fits the direction of our 1.8 to 3.5 fold inflation. It is a direction and not an explanation, and we claim no mechanism for the size.

The more useful result may be that the dish does well on a task closer to what it is used for. A drug-level sort into risk classes does not need FPD in milliseconds, and on that task the dish at its Cmax is about as good as the FDA's model at identifying low-risk drugs. The four-drug bridge, by contrast, asks the dish for a quantity it was never designed to give. We do not think the two results conflict. They say the assay is better as a ranking tool than as a measuring tool.

On individual prediction the evidence is weaker still. In the one study that grew cells from the same volunteers who took the drugs, subject-level dish slopes did not correlate with subject-level clinical slopes (Spearman ρ = −0.11, p = 0.70 for dofetilide). The informative control is that the same subjects' clinical responses to the two drugs did not correlate with each other either (ρ = 0.10, p = 0.73), so part of the weakness lies on the clinical side. Earlier studies reported positive correlations in other settings (Stillitano 2017; Shinozawa 2017). We use the bridge at the drug level only.

What is new here is modest, and mostly about method. We found no open implementation that runs the same correction cascade on dish and clinic data, carries provenance, and refuses to turn FPD into QT. We also found no public side-by-side at free-Cmax multiples that re-derives the clinical side from open ECG data, though effect-at-Cmax in the dish is not new. The quantification of informative missingness (1,887 of 2,067 nulls are EAD wells) is a measurement and not a discovery. The J-Tpeak exponent is not novel. These verdicts rest on a limited literature screen, and earlier screens in this project were overturned by sources the first pass missed, so none of them is a claim of priority.

## 7. Limits

The bridge rests on four drugs, so every statement about its size has wide uncertainty and none supports a calibration. The clinical population is healthy adults aged 19 to 35 with no arrhythmic events, which is the group least likely to show the outcome the dish is supposed to anticipate. The dish concentration axis is nominal, and this is the weakest link in the comparison. What is pipetted into a well is not what the ion channels see. Drug adsorbs to the plate plastic and to the perfusion tubing, and lipophilic compounds partition into the cell membrane and into intracellular lipid. This loss is separate from binding to serum protein in the medium, and it is compound-specific and saturable rather than a constant factor, so it does not cancel across a panel. Because the binding sites are finite, the fractional loss is largest at the lowest concentrations, which is where a careful assay works (Alvarez Baron 2025). The patch-clamp literature calls this nonspecific binding or drug loss rather than adsorption.

The size can be large enough to change a conclusion, and it varies widely by compound. In a study that verified concentrations by LC-MS/MS in the recording rig instead of trusting the nominal value, loss during hERG experiments ran from 6 % for ranolazine to 58 % for vanoxerine, with verapamil at 41 % and chloroquine at 19 %. Accounting for it moved the hERG IC50 to 2.2 nM for vanoxerine and to 0.75 uM for chloroquine (Alvarez Baron 2025). The same paper notes that earlier studies all plotted nominal concentrations, because none had measured the in-vitro concentration at all.

Two of those four compounds, ranolazine and verapamil, are bridge drugs here, and their losses differ sevenfold. That is the part that should worry a reader of our table: a per-compound loss does not cancel when two drugs are compared on a shared exposure axis. We cannot carry the numbers across, because those are patch-clamp rigs with flowing perfusion and ours are static MEA wells, so the mechanism mix differs: tubing and perfusion lines there, plate plastic and cell lipid here. What carries across is that the loss is real, compound-specific, and large enough to move a potency estimate.

Taken on its own this mechanism has a clear direction. If free concentration in a well is below nominal, every dish curve belongs further to the left on the exposure axis than we plotted it, so the cells are more potent per free molecule than they appear, and the dish over-read against the clinic at matched exposure would be larger rather than smaller. That is the direction of one mechanism in isolation, not a net correction, and we are careful here for a reason: a related claim about protein binding was made and withdrawn earlier in this project (section 5, item 1), because the one experiment that bears on the net sign ran the other way, with 8 of 18 drugs becoming arrhythmic at lower nominal concentrations when serum was present (Schocken 2018). Nobody measured free concentration in these wells, so the magnitude and even the net sign of the full correction are open. We are not reinstating the withdrawn claim. The in-silico comparison is out of sample for the 16 held-out drugs but qNet was developed with the 12 training drugs in view, and with n = 16 the comparisons have wide intervals. The coupled one-dimensional cable did not run: the OpenCL runtime on our remote machine exposes double precision on the device but not on the platform, myokit therefore used single precision, and the ORd model returned NaN. A single-cell model file also lacks the diffusion current binding that `set_conductance` needs, and in that state the call silently does nothing.

Several fits carry small technical limits too. The linear model under-covers (0.865 against 0.90), the degrees-of-freedom rule was not compared against Kenward–Roger, and a linear model on a curved verapamil response caused both misses against Vicente 2019. Residual diagnostics were not run.

## 8. What would settle the open questions

1. **Endpoint, site and medium.** A study in which the same cells and drugs are read by both electrodes and optical probes, in the same medium, at the same site.
2. **Hysteresis.** Beat-level or longer ECG recordings with the RR history, which public tables do not carry.
3. **The 6 to 7 ms shift.** Refit both studies with the full white-paper model (treatment indicator, nominal-time effects, baseline covariate), refit ECGCIPA day by day, and compare extraction methods, ideally on the waveforms of both.
4. **Missing wells.** A tipping-point analysis that assigns the missing wells a duration of the observed maximum plus δ and reports the δ at which a conclusion changes. This is not done.
5. **A usable bridge.** More single-agent clinical concentration–QTc data that are open. None exist publicly beyond these four.
6. **The cable model.** A computer with an OpenCL driver that exposes double precision at the platform level.
7. **The dish exposure axis.** Measurement of free compound in the well at the time of recording, by LC-MS/MS on the medium, for each drug and each site. Until that exists the dish axis stays nominal, and the size of the dish-clinic gap cannot be read as a property of the cells alone.

## 9. What this work is, and is not

None of us grew cells, recorded a field potential or ran a clinical study. All dish data come from the CiPA myocyte study and all clinical data from PhysioNet. The contribution is the curation and modelling half: assembling heterogeneous public data into an auditable form, implementing the correction cascade with provenance, fitting mixed-effects models that carry site and protocol variability explicitly, validating each stage against independent published analyses, decoding the clinical waveforms, running the in-silico arm and testing our own claims. Where the data could not answer a question, the product is a precise statement of why.

## 10. Challenges we met, and what they cost

These are the practical obstacles, as distinct from the scientific limits in section 7. Each one shaped a result, so each is worth stating.

**The source data arrive already corrected.** The published workbook ships ΔΔFPDc, not raw field potential durations and beat periods. Baseline, vehicle and rate correction were applied by the study, so our cascade cannot be exercised on the dish data. We verified the numbers instead (r = 0.999993 against all 420 published panels) and ran the same code path on an independent dataset (Millard 2018, 17 of 18 within 5 ms). The consequence is that the in-vitro baseline and vehicle steps are tested on synthetic data and on the clinical corpus, never on real raw wells.

**Two endpoints live in one column.** `ddFPDc` holds electrode FPD for 13 datasets and optical APD90 for 2. The paper calls them equivalent measures. As a scale they are not: the optical arm reads 1.44 times the consensus of the other nine sites. But the optical platform was used at one site, which is also the only serum-free site, and medium alone shifts baseline FPDc by 104 ms (Schocken 2018). Endpoint, site and medium move together, and no public dataset separates them.

**Missingness is informative and cannot be repaired.** 2,067 of 8,859 values are null and 1,887 of those are exactly the wells flagged as arrhythmic. The loss rises from 0.9 % to 60.8 % across the four concentrations. The underlying durations do not exist in the published data, so we bound the problem by fitting everything twice rather than imputing. A tipping-point analysis is the obvious next step and is not done.

**Small defects in the source tables.** Fifteen wells at the top concentration are under-coded for diltiazem, nitrendipine and verapamil, all at one site on one platform. One cell (clozapine, site 1, iCell², top concentration) has five measurable wells but no plotted value. Neither changes a conclusion; both are recorded rather than quietly dropped.

**One validation that still fails.** Study 9 (iCell² on CardioECR) is 22.0 ms off in the rate-correction check, with a beat-period coefficient of variation of only 0.061, so averaging order does not explain it. We could not account for it, and it is asserted in a test so it cannot be forgotten.

**A quality flag that does not flag the failure.** In the in-silico arm, beyond about 3.8 times free Cmax the model stops repolarising: the action potential duration is either missing or the detector returns a value measured from a depolarised resting potential. The run's own `ok` flag is true for all 408 runs, so the flag cannot be trusted and the physiology has to be checked directly. This was found while drawing Figure 2.

**Three technical blockers.** The coupled one-dimensional cable did not run: the OpenCL driver exposes double precision on the device but not on the platform, so the solver fell back to single precision and the model returned NaN. A single-cell model file also lacks the diffusion-current binding that the conductance setter needs, and in that state the call silently does nothing. Two papers behind paywalls could not be read, and one number in this report rests on an excerpt rather than a full text (Redfern 2003).

**Paper-to-paper disagreements that are not ours.** One published table disagrees with its own study's primary paper on the dofetilide effect; our value matches the primary paper to 0.2 ms. The two reference tables for free Cmax disagree for terfenadine by a factor of 13, which we traced to one of them being an exposure measured with a P450 inhibitor present.

**Our own errors.** Four claims were withdrawn after testing (section 5), and a citation in an earlier draft of this report carried numbers that had not been checked against the source. One of the three was wrong. It was corrected after reading the full text, and the episode is the reason every number here is read back from a file rather than quoted from memory.

## 11. Questions for the study authors

Several of the limits above are not limits of the public data so much as limits of what the data release records. The group that generated the CiPA myocyte study could answer most of them directly. The questions are ordered by how much they would change what can be concluded.

**On the pooled endpoint**

1. Did any site run both the electrode and the optical platform on the same cells, in the same medium? Even a small paired set would separate endpoint from site and medium, which no public dataset currently does.
2. Was the 1.44-fold higher reading at the optical site noticed during the study, and is it understood as an endpoint effect, a medium effect, or a site effect?
3. On what basis were FPD and APD90 treated as equivalent measures: agreement in rank, in drug effect, or in absolute scale? Our paired comparison gives optical = 1.51 × electrode + 24.8 ms.

**On the arrested wells**

1. For a well that develops an arrhythmic event, is any duration recorded before the event, or is the measurement abandoned at that point?
2. Is there a well-level record of why each value is missing, distinguishing arrhythmia from signal loss, from a well that never beat, from exclusion by a quality rule?
3. How were the missing wells handled in the published concentration-response analysis?

**On exposure**

1. Was the concentration in the well ever verified analytically, for any drug or any site? Recent patch-clamp work found losses from 6 % to 58 % between nominal and measured concentration.
2. Was anything done to limit adsorption to the plate and the tubing, and did the protocol differ between the serum-free site and the others?
3. Why does the free Cmax used in this study differ from the CiPA table for terfenadine by a factor of 13, and was the choice deliberate?

**On the protocol**

1. How long after drug addition was each reading taken, and did that interval differ by site? A published study reports beat rate rising 22.6 % and FPD falling 7.7 % over a 20-minute equilibration, which is the same order as the drug effects being measured.
2. What exactly do the baseline and vehicle corrections subtract: the same well before dosing, the plate mean, or the site mean, and over what time window?
3. Are the raw, uncorrected FPD and beat-period values retained anywhere? They would let the correction cascade be validated on real wells rather than on synthetic data.

**On specific results**

1. Does the shortening of FPD by verapamil appear at every site, and was it expected from its calcium-channel block? It is the one drug where the dish and the clinic disagree in sign, and a human ventricular model given verapamil's own channel profile agrees with the clinic.
2. Is anything known about study 9 (iCell² on CardioECR) that would explain a 22 ms discrepancy in a rate-correction check that the other 17 studies pass?
3. After the subject-matched result in 2019, which found no within-person correlation, is a larger subject-matched study planned? Our reading is that the clinical side contributed to that null, since the same subjects' responses to two drugs did not correlate with each other either.

## 12. The biology underneath: what the blockers do, and how the signals arise

This section is background. A reader who wants the physiology before the methods should start here.

### 12.1 The action potential, and why the dish's is different

A ventricular cell at rest sits near −85 mV, held there by the inward rectifier potassium current IK1. A stimulus opens fast sodium channels; sodium floods in and the cell depolarises within a millisecond or two (phase 0). Calcium then enters through L-type channels (ICaL) and holds the cell depolarised in a plateau, while potassium currents push the other way. The plateau ends when the potassium currents win: the rapid delayed rectifier IKr, carried by the hERG channel, and the slow one IKs. The cell returns to rest (phase 3), and IK1 holds it there (phase 4). APD90 is the time from the upstroke until the cell has recovered 90 % of the way back.

hiPSC-CMs are immature. They have little IK1, so they rest depolarised near −72 mV rather than −85 mV, and they beat spontaneously instead of waiting to be stimulated. That is not only a defect: it is a prerequisite for the assay, because a quiescent adult myocyte would have to be paced. But it changes the balance of currents, and that is the root of most of what follows.

### 12.2 What each blocker does

Using the measured block at each drug's free Cmax (`results/channel_block_at_crumb_free_cmax.csv`):

| Drug | hERG (IKr) | ICaL | Late sodium | Net effect |
|---|---|---|---|---|
| Dofetilide | 24.0 % | 0.4 % | 3.5 % | Almost pure IKr block: unopposed prolongation |
| Quinidine | 46.7 % | 8.1 % | 3.8 % | IKr block dominant, with Ito 13.9 % and IKs 8.3 % |
| Ranolazine | 22.6 % | — | 21.1 % | IKr block partly offset by late-sodium block |
| Verapamil | 14.6 % | 16.2 % | 0.5 % | IKr block offset by comparable calcium block |

Blocking IKr removes repolarising current, so the plateau lasts longer and the cell takes longer to recover. Blocking ICaL removes depolarising current, which shortens the plateau. Blocking late sodium does the same thing by a different route. A drug that does both — verapamil, ranolazine, mexiletine — has two effects that partly cancel, and which one wins depends on the cell.

This is why verapamil is the instructive failure. In an adult ventricular myocyte the two roughly balance and QTc barely moves (+3.5 ms in volunteers). In an immature cell with little IK1, the calcium current carries proportionally more of the plateau, so removing it shortens the field potential sharply (−86 ms in the dish). Running verapamil's own channel profile through a human ventricular model prolongs slightly, like the clinic — so the inversion is a property of the preparation, not of the compound.

### 12.3 From one cell to a field potential

A multielectrode array does not record the transmembrane potential. The electrodes sit outside the cells, under a sheet of them, and record the extracellular field that the sheet generates — closer to the time derivative of the membrane potential than to the potential itself. The sharp downward spike is the synchronised sodium upstroke. The slow hump much later is the repolarisation wave, and its peak is the conventional end of the field potential duration. FPD therefore tracks repolarisation and is sensitive to IKr block, but it is not an action potential duration and the two are not interchangeable (section 4.4 and `results/bland_altman_stats.csv`).

### 12.4 From many cells to an ECG

The ECG is a whole-heart signal and its shape comes from *differences* between cells, not from any one cell. The QRS complex is the wave of depolarisation sweeping through the ventricles, carried by the His-Purkinje system. The T wave is the spread of repolarisation times across the wall and from apex to base: if every cell repolarised simultaneously the T wave would vanish. Mid-myocardial (M) cells have the longest action potentials, so the end of the T wave is set largely by the last tissue to recover.

The QT interval is therefore activation time plus repolarisation time of the latest region. Two consequences matter here. First, electrotonic coupling between neighbouring cells blunts the differences: real tissue shows about 30 to 40 ms of dispersion where uncoupled cells would show far more. In our uncoupled simulations the QT surrogate is simply the M-cell APD90 in 121 of 136 conditions, with an implied dispersion of 101.9 ms — which is why the common "weighted average APD90" mapping has no mechanism behind it. Second, a dish monolayer has no activation sequence, no Purkinje system, no transmural gradient and no autonomic input, so it has no analogue of the T wave at all.

## 13. The four corrections, and what each one is doing

### 13.1 Rate correction

Repolarisation takes longer when the heart beats slower, so a raw QT or FPD is uninterpretable without knowing the rate. The standard repair divides by a power of the cycle length: Ic = I / CL^α, with α = 1/3 for Fridericia and 1/2 for Bazett. We fit α on drug-free records only — baseline and placebo — because fitting on treated records would absorb part of the drug effect into the correction. The fitted exponents came out at 0.301 and 0.299, close to Fridericia's 1/3, which is both the expected answer in healthy volunteers and a check on the code.

The dish needs this more than the clinic does, not less: hiPSC-CMs beat spontaneously at roughly 29 to 37 beats per minute, and a drug that changes the beat rate changes the raw FPD for reasons that have nothing to do with repolarisation.

### 13.2 Hysteresis

QT does not follow a change in heart rate immediately. It lags, with a time constant of roughly one to two minutes, so the QT in any given beat reflects not that beat's RR but a weighted average of the preceding minute or two. The standard correction forms exactly that: an exponentially weighted average of preceding RR intervals, used in place of the instantaneous RR.

**We could not do it, and the reason is structural rather than a matter of effort.** Each public record is a 10-second strip. That is 8 to 17 % of one adaptation time constant, so most of the RR history that determines the QT in the strip happened before the strip began. No weighting of 10 seconds of data can recover information that is not in it. The cost is measurable: within a single window, QTcF computed on the slowest and the fastest beat differs by a median of 12.6 ms, and by more than 10 ms in 66.8 % of 238 records. The quantity being averaged moves by more than the regulatory threshold inside one recording.

### 13.3 Autonomic tone

There is no standard "autonomic correction", and it is worth being clear about why. Sympathetic and parasympathetic traffic change both the heart rate and the repolarisation itself, and the second part is not removed by rate correction. In practice the field controls it by design rather than correcting it: recordings supine, at matched times of day, fasted, with a placebo arm subtracted so that circadian and postural common-mode effects cancel.

The size of what is being controlled is not small. Giving ibutilide to healthy subjects with intact autonomic tone prolonged QTc by 43 ± 8 ms; under double autonomic blockade the same drug produced a further 68 ± 5 ms. Tone *attenuates* drug-induced prolongation. A dish has no autonomic input at all, so it sits permanently in the state where a hERG blocker's effect is largest. That is a structural difference between the systems, not something either side can correct away, and it runs in the right direction to explain part of the dish's over-reading.

### 13.4 Irregular beats, and the median beat

QT is only defined on a sinus beat that follows a reasonably regular preceding beat. An ectopic beat, or the compensatory pause after it, changes repolarisation for reasons unrelated to the drug. The usual handling is to exclude non-sinus beats and require a stable preceding RR.

The **median beat** is the standard device: align every beat in the window on a fiducial point, usually the R peak, then take the median across beats sample by sample. It suppresses noise and is robust to a single aberrant beat, which is why it is preferred to the mean.

One correction to a claim that circulates: the median beat does **not** deliver the √N improvement in signal-to-noise that averaging N beats would. That property belongs to the mean. The median is roughly 1.25 times less efficient than the mean for Gaussian noise. It is used for robustness, not for efficiency.

In the dish the analogous problem is more severe and is not solvable by averaging. A well that develops an early afterdepolarisation has no measurable field potential duration at all, so the observation is lost rather than noisy — and it is lost precisely when the drug effect is largest. That is the informative missingness that rises from 0.9 % to 60.8 % across the concentration range.

## 14. Why one solution does not fit all three arms

The cascade is applied symmetrically, which is the point of the package. But symmetry of *form* is not symmetry of *difficulty*. Every stage means something different in each arm, and several problems admit a fix in one arm that is meaningless or impossible in another. This is the part most easily lost when the three are described as parallel.

| Stage | In vitro | In vivo | In silico |
|---|---|---|---|
| **Rate** | The preparation sets its own rate, 29–37 bpm, and drugs change it. Correction is large and unavoidable. | Rate is set by autonomic tone; correction is modest at resting rates. | Rate is imposed by the pacing protocol. There is nothing to correct. |
| **Baseline** | The well before dosing. Available in principle, but **not shipped** in the public workbook. | The subject's own pre-dose record, within period. Available. | The drug-free run of the same model. Exact, and noise-free. |
| **Control** | Vehicle wells, usually DMSO. | A placebo arm, time-matched. | None needed, and none exists. |
| **Hysteresis** | Not applicable: the preparation is at steady state at its own rate. | Applicable and important, but **unidentifiable** in 10-second records. | Not applicable: pre-paced to steady state by construction. |
| **Autonomic tone** | Absent entirely. Cannot be corrected, only acknowledged. | Present, controlled by design rather than by correction. | Absent by construction, which is a modelling choice rather than a limitation. |
| **Missing data** | Informative and severe: arrested wells, 0.9 % → 60.8 % with concentration. | Rare and close to at-random. | Informative but deterministic: beyond about 3.8× the model stops repolarising. |
| **Exposure** | Nominal concentration, reduced by protein binding, plastic and membrane adsorption, by an unmeasured amount. | Measured free plasma concentration. | Exactly the concentration prescribed. No uncertainty at all. |
| **Dominant variance** | Between site × platform × cell line. | Between subjects. | Parameter uncertainty — and only for hERG. |
| **Endpoint** | FPD, or APD90 in the optical arm, pooled in one column. | QT, J-Tpeak, Tpeak–Tend. | APD90, qNet. |

Four consequences follow, and each one blocks a repair that works elsewhere.

**Averaging fixes noise in the clinic and loses data in the dish.** A noisy clinical beat can be averaged away with a median beat. An arrhythmic well produces no number at all, so there is nothing to average; the observation is destroyed exactly where the drug effect is largest.

**Placebo subtraction removes common-mode drift in the clinic and cannot exist in silico.** The simulation has no circadian rhythm, no meal, no posture change, so the correction that rescues the clinical estimate has no referent. Applying it would be meaningless rather than conservative.

**Hysteresis is the clinic's problem alone.** The dish and the model are at steady state. Spending effort on an effect-compartment model would improve neither.

**Exposure uncertainty runs in opposite directions.** In silico the concentration is known exactly and the *response* is uncertain. In the clinic the concentration is measured with error but the response is observed. In the dish the concentration is the least certain quantity in the whole chain, and it is the one quantity on which the entire comparison rests. That asymmetry is why the shared axis is the weakest link, and why it gets its own paragraph in section 7.

The practical conclusion: only the **concentration axis** is genuinely shareable. Everything else must be handled per arm and then compared, which is what the cascade does, and it is why a single "harmonisation" procedure applied to all three would be wrong even if it were convenient.

## 15. Is "are the three systems consistent?" a sensible question?

Asked plainly, no. Asked precisely, yes — and the precise version is the one the field has standardised on.

**Why the plain version fails.** The three measure different observables, so agreement in milliseconds is undefined. Consistency also has no value without an exposure attached: nonclinical assays show robust specificity but poor sensitivity for clinical QTc prolongation at 1–30× clinical exposure multiples, with sensitivity rising to 0.67–0.80 at higher multiples for hERG and in vivo QTc. The same assay is therefore concordant or discordant depending only on where it is evaluated. And the reference standard is itself a surrogate: a thorough-QT outcome measured in healthy volunteers who experience no arrhythmia.

**What the field measures instead.** Diagnostic accuracy against a named clinical endpoint — sensitivity, specificity, predictive values, likelihood ratios — rather than agreement between assays. The most useful form is Bayesian: against an estimated 20 % pretest probability, a positive hERG result raises the post-test probability of clinical QTc prolongation to 53 %, a negative result lowers it to 9 %, and a double-negative nonclinical result lowers it to 4 %. The ICH E14/S7B Q&A builds on exactly that asymmetry, using a negative hERG assay together with a negative in vivo study to support the absence of a clinically relevant QT liability.

Note what the regulation rests on: two *negatives* agreeing. Not on the systems producing matching numbers.

**Would a relative scale rescue the plain version?** Partly, and the partial failure is informative. We tested it directly (`results/relative_scale_test.csv`). On a percent-of-baseline scale the dish-to-clinic gap narrows from 1.81× to 1.48× but does not close. Normalised instead to a reference drug measured in the same system, the two *pure* hERG blockers agree almost exactly — dofetilide-normalised quinidine is 1.023 in the dish against 1.025 in the clinic, a ratio of 0.999 — while ranolazine is 1.92× off and verapamil inverts. So a relative scale removes the units problem and the common multiplicative factor, and leaves the mechanism-specific disagreement untouched. The residue is not a scaling artefact and no normalisation will remove it.

Two cautions on that test: the dish baseline a percent scale requires is not in the public workbook and had to be imported from an independent study, and verapamil's clinical interval spans zero, so any ratio with it in the denominator has unbounded variance.

**And would a Bland-Altman plot settle it?** Only where the two methods measure the same quantity. For dish against clinic it is not licensed: the average of an FPD and a QT estimates nothing, and four drugs cannot support limits of agreement. Where it *is* licensed, it is informative (`results/bland_altman.png`). Optical APD90 against electrode FPD, over 97 paired drug × concentration cells, gives a bias of +34.7 ms and limits of agreement from −227 to +297 ms, with strong proportional bias (slope 0.800, p < 0.001; 0.502 with the single extreme cell removed). That 524 ms-wide band, against drug effects of 10 to 130 ms, is the sharpest available test of the source study's description of the two as equivalent measures. Our tangent T-end against the study's own, on 220 records in 22 subjects, gives a bias of −9.1 ms with limits of −64 to +46 ms and no proportional bias.

Limits of agreement there use the closed-form unbalanced one-way variance components rather than a naive SD, because several concentrations share a drug and several records share a subject; treating them as independent would have understated the spread. A random-intercept model was tried first and did not converge on the clustered comparison.

**The useful form of the question** is therefore not "are they consistent" but "what does each one add that the others do not". Three systems that always agreed would make two of them redundant. The framework pays for the combination precisely because they fail in different directions.

## 16. What the data cannot support, and what we would want

### 16.1 Random effects that are not identifiable

**Below the dataset level, in vitro.** The workbook gives one ΔΔFPDc per well and no plate identifier, well position, electrode count, timestamp, cell lot or passage number. Plate-to-plate variation and edge effects — both well documented in plate assays — therefore have nowhere to go in the model. The finest grouping the data supports is dataset = site × cell line × platform, which is what we fit; everything below it falls into residual. This matters for interpretation: the dish intervals are 25 to 106 ms wide and we attribute that to between-dataset heterogeneity, but we cannot exclude plate structure as part of it.

**Across channels, in silico.** `data/hill28.csv` carries IC50 and Hill coefficient as single point estimates for all seven currents. Only hERG has a distribution (`data/herg_boot28.csv`, 20 bootstrap sets per drug drawn from the 2,000 published). Every qNet interval in this work therefore propagates hERG binding uncertainty alone and treats the other six channels as known exactly. The true uncertainty on qNet is larger than reported, by an amount the public data cannot reveal.

### 16.2 The data we would ask for, in order of value

1. Raw, uncorrected well-level FPD and beat period, with plate and well coordinates and read timestamps. This one item unlocks cascade validation on real wells, the plate random effect, percent-change normalisation, and the equilibration-drift question at once.
2. Measured free concentration per well, per drug, per site.
3. A paired endpoint set: the same cells read by both platforms in the same medium.
4. A reason code for every missing well, separating arrhythmia from signal loss, from never-beating, from quality exclusion.
5. IC50 and Hill uncertainty for all seven channels, not hERG alone.
6. Clinical recordings longer than 10 seconds, carrying the preceding RR history, so hysteresis becomes identifiable.
7. Beat-level clinical annotations rather than one median beat per record.
8. More open single-agent concentration–QTc datasets. Four drugs is the entire bridge.
9. A subject-matched dish-and-clinic study larger than 16 people.
10. Any torsades outcomes at all. The corpus contains none.

## 17. Literature summary, paper by paper

What each key paper contributes, what we used it for, and the one caveat that matters when leaning on it. Papers marked **(not read in full)** were used only via abstract or search result and are flagged as such wherever they are cited.

### The in-vitro evidence

**Blinova 2018, *Cell Reports* 24:3582–3592** — the CiPA myocyte study. Twenty-eight drugs across ten sites, two commercial hiPSC-CM lines and six platforms: the largest public dataset of its kind, and the entire dish side of this work. *Caveat:* it ships ΔΔFPDc only, already baseline-, vehicle- and rate-corrected, and it pools electrode FPD with optical APD90 in that one column while describing them as equivalent measures. Section 4.4 and our Bland-Altman analysis both bear on that description.

**Blinova 2016, *Toxicological Sciences* 155:234–247** — the earlier pilot, which established that the assay separates torsadogenic from non-torsadogenic compounds categorically. Useful as the precedent for drug-level use; not used for any number here.

**Blinova 2019, *Clinical and Translational Science* 12:687–697** — "clinical trial in a dish". Cells reprogrammed from 16 of the volunteers who took the drugs, so the dish and the clinic refer to the same people. Found no within-person correlation (ρ = −0.11, p = 0.70 for dofetilide). *The informative control:* those subjects' clinical responses to two drugs did not correlate with each other either (ρ = 0.10), so part of the null is on the clinical side. This is why the package stays at drug level.

**Millard 2018, *Toxicological Sciences* 164:550–562** — a blinded cross-site MEA pilot whose supplementary table reports baseline beat period, FPD and corrected FPD for 18 studies. The only public dataset on which our rate correction could be validated against an independent implementation: 17 of 18 within 5 ms, median 0.65 ms. *Caveat:* study 9 is off by 22.0 ms and we cannot explain it.

**Schocken 2018, *J Pharmacol Toxicol Methods* 90:39–47** — media effects. Serum-containing versus serum-free medium changes baseline FPDc by 104 ms (575 ± 59 against 471 ± 44) and slows beating by a third. *Why it matters here:* it is the reason the optical-site confound cannot even be signed. It also overturned one of our claims — of 18 drugs that induced arrhythmias, 8 did so at *lower* concentrations with serum present, the opposite of what protein binding alone predicts.

**Guerrelli 2024, *Am J Physiol Heart Circ Physiol* 327:H12–H27** — temporal drift. Beat rate rose 22.6 % and FPD fell 7.7 % over a 20-minute equilibration. The same order as the drug effects being measured, which is why read-out timing is one of our questions for the study authors.

**Ernault 2024, *Am J Physiol Heart Circ Physiol* 326:H800–H811** — what a field potential actually is, and why the criteria for reading it differ from those for an electrogram. The reference for section 12.3.

### The clinical evidence

**Johannesen 2014, *Clin Pharmacol Ther* 96:549–558** — the ECGRDVQ trial: dofetilide, quinidine, ranolazine, verapamil and placebo in 22 subjects, crossover, with concentrations. Three of the four bridge drugs come from here.

**Johannesen 2015, *Clin Pharmacol Ther* 99:214–223** — ECGDMMLD, late-sodium block. We reproduce its dofetilide-alone result to within 0.2 ms at two timepoints. *Why it matters:* it resolved what had looked like our own 6–7 ms error, which turned out to be a disagreement between two published analyses of the same trial.

**Vicente 2015, *J Am Heart Assoc* 4** — T-wave morphology in ECGRDVQ; its Table 2 is our primary validation target, 8 of 8.

**Vicente 2016, *PLoS ONE* 11:e0163619** and **Vicente 2019, *Clin Pharmacol Ther* 105:943–953** — the late-sodium biomarker work and the CiPA phase I validation study. Vicente 2019 is the source of the ECGCIPA data and of the prespecified J-Tpeak power-law exponent of 0.58, which we confirmed on 5,740 records. *Caveat:* its Table 3 and the primary paper disagree on dofetilide; we match the primary paper.

**Johannesen 2016, *PLoS ONE* 11:e0166925** — the automated J-Tpeak and Tpeak–Tend algorithm behind the interval tables we consume.

**Hnatkova 2019, *Scientific Reports* 9** — the linear J-Tpeak rate correction, the alternative to the power law. Per record the two differ by a median of 0.44 ms; at drug level by 1.2 ms or less.

**Malik 2019, *Drug Safety* 42:401–414** — why population rate corrections are inadequate and individual QT/RR profiles matter. The basis for the per-subject exponent option.

**Garnett 2016, *Clin Pharmacol Ther* 100:170–178** — operating characteristics of linear concentration–QT models, the design our clinical fits follow.

**Strodthoff 2023, *Scientific Data* 10** — PTB-XL+, where three algorithms measured the same 21,795 ECGs. Our external evidence that interval extraction carries a rate-dependent bias of about 6 ms across a 300 ms RR change.

**Smith 2007, *J Cardiovasc Electrophysiol* 18:960–964** — autonomic tone attenuates drug-induced QT prolongation. Ibutilide prolonged QTc by 43 ± 8 ms with tone intact and by a further 68 ± 5 ms under double blockade. The quantitative basis for expecting a denervated preparation to over-read.

### The in-silico and ion-channel evidence

**O'Hara 2011, *PLoS Comput Biol* 7:e1002061** — the human ventricular model used for our 408 simulations. The paper reports endocardial APD of about 275 ms in human preparations at a 1,000 ms cycle length, matched by the model at 273 ms; our control endocardial APD90 is 265.5 ms, 2.7 % below. *Caveat:* epicardial and M-cell control values are not stated numerically in the text, so those remain unverified against it.

**Dutta 2017, *Frontiers in Physiology* 8:616** — the ORd-CiPA reparameterisation and the origin of qNet. **Li 2017, *Circ Arrhythm Electrophysiol* 10** — hERG binding kinetics combined with multichannel pharmacology, the model behind the trapping formulation. **Li 2019, *Clin Pharmacol Ther* 105:466–475** — the formal assessment of the in-silico metric under CiPA.

**Crumb 2016, *J Pharmacol Toxicol Methods* 81:251–262** — 30 drugs against the CiPA channel panel. Table 1 is our source of free Cmax. *Caveat we had to correct:* it does **not** tabulate protein binding, contrary to an earlier note in this repository.

**Alvarez Baron 2025, *Clin Transl Sci* 18:e70451** — patch-clamp protocol practice, including concentration verification in the recording rig. Losses from nominal ranged from 6 % (ranolazine) to 58 % (vanoxerine), with verapamil at 41 % and chloroquine at 19 %. The quantitative anchor for section 7's treatment of the dish exposure axis. *Caveat:* perfusion rigs, not static MEA wells; the numbers do not transfer.

**Pritchard 1985, *J Clin Pharmacol* 25:347–353** — equilibrium-dialysis measurement of bepridil plasma binding, and the evidence that its binding saturates.

### The framework and translation literature

**Sager 2014, *Am Heart J* 167:292–300**; **Colatsky 2016, *J Pharmacol Toxicol Methods* 81:15–20**; **Fermini 2016, *SLAS Discov* 21:1–11** — the CiPA paradigm, its rationale and progress. Context rather than data.

**Park 2018, *British Journal of Pharmacology*** — the HESI consortium analysis of whether nonclinical repolarisation assays predict thorough-QT outcomes, using sensitivity, specificity, ROC curves, predictive values and likelihood ratios. The source of the finding that specificity is robust while sensitivity is poor at 1–30× exposure multiples. **(not read in full)**

**Gintant 2022, *British Journal of Pharmacology*** — the Bayesian framing: pretest 20 %, positive hERG 53 %, negative hERG 9 %, double-negative 4 %. The clearest statement of why "consistency" is the wrong target and incremental information is the right one. **(not read in full)**

**Strauss 2021, *Clin Pharmacol Ther*** — the ICH E14/S7B updates and the discriminatory power of the core assays for clinical QTc and for TdP. Notes that exact probabilities should be cited with caution. **(not read in full)**

**Stillitano 2017, *eLife* 6** and **Shinozawa 2017, *Stem Cell Reports* 8:226–234** — the two studies reporting positive subject-level correlations, against which Blinova 2019's null must be read. The literature is unsettled and this package does not try to settle it.

**Redfern 2003, *Cardiovasc Res* 58:32–45** — the provisional safety-margin framework. Closed access; one excerpt was read from a search result and is the basis for the terfenadine exposure note. **(not read in full)**

## 18. Were the processing choices right? A methodological audit

Every number here depends on four choices that could have gone another way: the rate-correction exponent, the shape of the concentration predictor, the random-effect structure, and whether the placebo adjustment is pre-formed or estimated inside the model. This section states each model in full, refits the headline quantity under the alternative, and reports what moves (`results/methodology_sensitivity.csv`, `scripts/methodology_sensitivity.py`).

### 18.1 The models, written out

**In vivo, as implemented here.** For subject *i* and post-dose record *j*, with concentration C centred at the geometric-mean observed Cmax C\*:

    ΔΔQTcF_ij = β₀ + β₁(C_ij − C*) + b_i + ε_ij,   b_i ~ N(0, σ²_b),  ε_ij ~ N(0, σ²)

β₀ is the effect at C\*, reported with a two-sided 90 % interval using a t quantile on (number of subjects − 1). The double delta is formed *before* the model: baseline subtracted within subject and period, then the time-matched placebo arm subtracted.

**In vivo, as the field prescribes.** The prespecified linear mixed-effects model of the concentration-QTc white paper uses the single delta as the response and estimates the placebo adjustment inside the model:

    ΔQTcF_ij = θ₀ + θ₁C_ij + θ₂(QTc0_i − mean QTc0) + θ₃TRT_ij + Σ_k θ_4k·1{time = k} + b_i + b'_i C_ij + ε_ij

with placebo records entering at C = 0, TRT = 0, and the reported quantity ΔΔQTc(C) = θ₁C + θ₃. In that model, study treatment and post-baseline time point are fixed effects with random effects on intercept and slope per subject, and placebo adjustment is achieved through the treatment and time-point effects. The centred baseline QTc enters as a covariate.

**In vitro.** For dataset *d* (site × cell line × platform, 15 of them) and well *w*, with m the multiple of free Cmax:

    ΔΔFPDc_dw = γ₀ + γ₁(log₁₀ m_dw − log₁₀ m*) + u_d + ε_dw,   u_d ~ N(0, τ²)

The random effect is the dataset, not the well, because nothing below the dataset level is identifiable (section 16.1).

**In silico.** No mixed model, and this is not an omission. The simulation has no sampling error: a given parameter vector returns one answer. What uncertainty exists is parametric, y_b = g(θ_b) over bootstrap draws θ_b of the hERG binding parameters, and it covers one channel of seven. The band in Figure 2c is the range across the three cell types — a spread, not a confidence interval — and the figure says so.

### 18.2 Does one cascade make sense?

In form, yes; in parameters, no. The same ordered sequence — rate, baseline, control, model, shared scale — is what makes the comparison auditable, and enforcing it in code is what stops an endpoint being used before it has been corrected. But every parameter inside it has to be arm-specific, and two of the five steps are vacuous in at least one arm (section 14). A single cascade with a single set of parameters would be wrong.

### 18.3 What if the power law differs by arm?

It does, and this is where the largest mistake could have been made.

**Clinically**, Fridericia is defensible on our own data: exponents fitted on drug-free records come out at 0.301 and 0.299, close to 1/3. The sensitivity analysis shows why the choice matters anyway:

| Drug | Fridericia (used) | Bazett | Individual exponent |
|---|---|---|---|
| Dofetilide | 73.1 [68.5, 77.7] | 72.4 [67.3, 77.5] | 72.3 [67.4, 77.3] |
| Quinidine | 74.9 [69.3, 80.5] | 85.4 [78.4, 92.4] | 78.6 [72.1, 85.1] |
| Ranolazine | 11.2 [7.6, 14.8] | 14.4 [9.6, 19.2] | 12.4 [7.8, 17.0] |
| **Verapamil** | **3.5 [−0.1, 7.0]** | **10.5 [6.0, 15.0]** | 5.8 [1.7, 9.8] |

**Bazett flips the ICH E14 decision for verapamil.** Under Fridericia the upper bound is 7.0 ms and the drug is negative; under Bazett it is 15.0 ms and the drug is positive. This is the textbook failure of Bazett — it over-corrects when heart rate moves, and verapamil slows the heart — but it is worth seeing it change a regulatory call on real data rather than reading about it. The result is asserted in the test suite so it cannot be quietly lost.

**In vitro** we could not make this choice, because the workbook arrives already corrected. Fridericia is the field standard for field potentials: the Axion induced-model standard defines FPDc = FPD/Beat Period^0.33, and the convention is adopted to minimise the influence of a wide range of interspike intervals. But the standard is not unchallenged. One characterisation study notes that there is no evidence that clinical correction formulae are applicable to FPD from hiPSC-CMs of varying beat rate, finds both Fridericia and Bazett fairly effective, and then reports that Fridericia likely under-corrects and Bazett over-corrects at slow interspike intervals. A dedicated optimisation study developed and qualified correction algorithms using the If inhibitor ZD7288 and isoproterenol to move the beat rate deliberately, and found that after Fridericia's correction the number of plates with significantly different FPD values was reduced.

So the best available method in vitro is **not** to assume an exponent but to qualify one per preparation with a rate-changing tool. We could not do that. The only check available was indirect: the exponent implied by the Millard 2018 baseline data is 0.324 against Fridericia's 0.333, which is reassuring but is a different study's cells.

### 18.4 What if the predictor shape were the same in both arms?

This is the choice with the largest consequence, and it cannot be shared.

| Drug | log₁₀ multiples (used) | linear in multiples |
|---|---|---|
| Dofetilide | +132.0 [93.1, 170.8] | +96.2 [60.9, 131.6] |
| Quinidine | +121.7 [76.6, 166.7] | +126.8 [85.5, 168.1] |
| Ranolazine | +39.5 [26.4, 52.6] | +33.8 [21.2, 46.4] |
| Verapamil | −87.2 [−99.3, −75.1] | −43.4 [−53.8, −33.0] |

Forcing a linear predictor halves verapamil's effect and cuts dofetilide's by a quarter, because tested in-vitro concentrations span up to 3,500-fold and a linear fit is dominated by its top point. Conversely, forcing a log predictor on the clinical side would be wrong: clinical concentrations span less than two decades, and the ICH E14 decision is defined on a linear concentration scale. The two arms need different predictors for reasons internal to each, which is exactly what "one cascade" must not be taken to mean.

### 18.5 What if we pooled differently?

**Across datasets, in vitro.** Treating the dataset as a fixed rather than a random effect moves dofetilide from +132.0 to +118.3 ms, quinidine from +121.7 to +123.6, ranolazine from +39.5 to +39.3, and verapamil from −87.2 to −86.2. The random-effect version shrinks datasets toward the common mean and supports a statement about sites in general; the fixed-effect version describes these fifteen datasets and nothing else. We want to generalise, so the random effect is the right default — but the dofetilide gap of 13.7 ms shows the choice is not free.

**Across studies, clinically.** We deliberately do not pool ECGRDVQ with ECGCIPA, because at matched concentration they differ by 6 to 7 ms for the same drugs and we cannot explain it (section 4.7). Pooling would average that discrepancy into a number that looks more precise than the evidence warrants.

### 18.6 Where our clinical model deviates from the standard, and what it costs

Four deviations, honestly stated.

1. **We pre-form the double delta** instead of estimating the placebo adjustment inside the model through treatment and time effects. Tested directly: the white-paper-style model gives 69.8 against our 73.1 (dofetilide), 71.7 against 74.9 (quinidine), 11.7 against 11.2 (ranolazine) and 3.5 against 3.5 (verapamil). **The largest gap is 3.4 ms and no decision changes.**
2. **No centred baseline QTc covariate.** The white paper requires it and we omit it. This is a genuine gap, not a tested-and-rejected choice, and it is the deviation we would fix first.
3. **Random intercept only.** Adding a random slope moves the effects by 0.0 to 3.0 ms and changes no decision; a separate earlier test over nine published comparisons moved the mean absolute gap from 3.16 to 3.04 ms.
4. **No categorical time effect.** Partly defensible, since time-matched baseline adjustment is argued to minimise diurnal variation so that categorical time effects are not needed, but we did not test it.

**The best method, if starting again**, is the prespecified model as written: ΔQTcF as the response, concentration plus centred baseline plus treatment plus time as fixed effects, random intercept and slope per subject, and the double delta derived from model contrasts rather than formed in advance. We now know what the shortcut cost on these four drugs — at most 3.4 ms, with no change to any conclusion — which is the honest way to report a simplification: not by claiming it does not matter, but by measuring how much it does.


### 18.7 Is the concentration-response linear at all?

**The in-silico arm fits nothing.** This is worth saying plainly, because Figure 2c looks like a steeply curving dose-response and invites the question. The simulation returns one action potential per parameter vector; there is no sampling, no residual and no regression. Where we quote an in-silico effect at 1.07x free Cmax, that multiple was simulated directly, so even interpolation is not involved.

**And the curve in that figure is partly a log-axis artefact.** Plotted against concentration rather than log concentration, the simulated response is close to a straight line across the exposures that matter: a linear fit of dAPD90 on multiples of free Cmax gives R^2 = 0.9993 for dofetilide (largest residual 2.2 ms), 0.9894 for quinidine, 0.9820 for ranolazine and 0.9089 for verapamil. The convexity is real but it is a supratherapeutic phenomenon: extending dofetilide to 30x in the endocardial cell, which does keep repolarising, the linear fit falls to R^2 = 0.9387 with a largest residual of 70 ms. So the model is near-linear through the therapeutic range and accelerates well above it, which is also where it eventually fails to repolarise.

**Where linearity *is* assumed, we tested it.** The clinical model is linear in concentration. That is the field default, and the same guidance that prescribes it requires the assumption to be checked and an Emax model fitted when it fails; study protocols routinely prespecify exactly that fallback. Fitting linear, log-linear and Emax forms with the same random intercept, by maximum likelihood so the fixed structures are comparable (`results/model_shape.csv`):

| Drug | Best by BIC | linear ΔBIC | effect at Cmax, linear | effect at Cmax, best | E14 call changes? |
|---|---|---|---|---|---|
| Dofetilide | **linear** | 0.0 | +73.1 [68.6, 77.6] | +73.1 | no |
| Quinidine | log-linear | 19.3 | +74.9 [69.4, 80.4] | +73.2 [68.0, 78.4] | no |
| Ranolazine | log-linear (tie) | 0.4 | +11.2 [7.7, 14.8] | +11.6 [8.1, 15.2] | no |
| Verapamil | Emax | 7.1 | +3.5 [0.0, 6.9] | +4.0 [0.6, 7.4] | no |

So linear is the best form for only one of the four. Quinidine clearly prefers a curved form, which is independently corroborated: a recent re-analysis of quinidine concentration-QT data found the Emax and log-linear models were better fits than the linear model by BIC. Verapamil prefers Emax, which is also the drug whose two misses against the published analysis we had already attributed to fitting a linear model to a curved response.

**But the consequence is small.** The effect at the geometric-mean Cmax moves by at most 6 ms across all three shapes, and **no ICH E14 decision changes for any of the four drugs**. Compare that with the rate-correction choice in section 18.3, where swapping Fridericia for Bazett flips verapamil's decision outright. Curvature is the more interesting modelling question and the less consequential one.

**The alternatives, in the order we would try them.** Emax, f(C) = Emax·C/(C + EC50), is the standard first fallback and is what protocols prespecify. Sigmoid Emax adds a Hill coefficient and is the form the underlying receptor pharmacology actually implies. Log-linear, f(C) = p1·log(1 + C/p2), is a convenient middle ground. For the in-vitro arm we already use a nonlinear transform — the predictor is log10 of the multiple — which is a local approximation to the middle of a sigmoid over a 3,500-fold concentration span.

**The larger gain is elsewhere.** The same re-analysis found that adding hysteresis to the linear model improved BIC by about 110, against the 4-point gap between Emax and linear. The shape of the concentration axis is a second-order question next to the lag between concentration and effect — and hysteresis is exactly what the 10-second public records cannot support (section 13.2). If one modelling improvement were available to this project, it would not be Emax. It would be longer recordings.


## 19. Availability

The code is the `dish2ecg` package (133 tests, version 0.2.2). `results/README.md` maps every result file to the script that produced it and gives the regeneration order. The repository is published at https://github.com/rejaul/dish2ecg. It is not pushed from this environment, which has no GitHub credential configured; the archive carries its own git history so a single `git push` completes it. The well-level dish table is not redistributed; `scripts/build_wells.py` rebuilds it from the original workbook.

## Glossary

- **QT interval.** The time on an ECG from the start of ventricular depolarisation to the end of the T wave. Long QT predisposes to torsades de pointes.
- **QTc, QTcF.** QT corrected for heart rate; the F is Fridericia's correction (divide by RR to the power 1/3).
- **FPD.** Field potential duration, the dish analogue measured with electrodes under the cells.
- **APD90.** Action potential duration to 90 % repolarisation, measured in single cells or with optical probes.
- **Δ and ΔΔ.** Change from baseline; change from baseline minus the placebo or vehicle change.
- **EAD.** Early afterdepolarisation, an abnormal extra depolarisation during repolarisation. In the dish it is a marker of arrhythmic behaviour.
- **Free Cmax.** The peak concentration of unbound drug in plasma at a clinical dose.
- **Nominal concentration.** In a dish experiment, the concentration calculated from the amount of drug added to the well. The concentration actually available to the cells is lower, because drug binds medium protein, sticks to plastic and tubing, and partitions into the cells themselves.
- **Fraction unbound (fu).** The share of drug in plasma that is not bound to protein.
- **CiPA.** Comprehensive in vitro Proarrhythmia Assay, the regulatory initiative that proposed combining ion-channel data, a computer model and hiPSC-CM data.
- **qNet.** The CiPA in-silico score: net charge carried by the main ion currents during a beat, in µC/µF. Lower values mean higher torsades risk.
- **Hysteresis.** The delay with which QT follows a change in heart rate.
- **AUC.** Area under the ROC curve; 0.5 is chance, 1.0 is perfect separation.

## References

- Alvarez Baron C, et al. (2025) Reassessing the conduct of patch clamp cardiac ion channel assays to improve nonclinical data translation to clinical ECG changes and proarrythmia risk. *Clin Transl Sci* 18:e70451. doi:10.1111/cts.70451
- Blinova K, et al. (2016) Comprehensive translational assessment of human-induced pluripotent stem cell derived cardiomyocytes for evaluating drug-induced arrhythmias. *Toxicol Sci* 155:234–247. doi:10.1093/toxsci/kfw200
- Ernault AC, et al. (2024) Interpretation of field and LEAP potentials recorded from cardiomyocyte monolayers. *Am J Physiol Heart Circ Physiol* 326:H800–H811. doi:10.1152/ajpheart.00463.2023
- Garnett C, et al. (2016) Operational characteristics of linear concentration-QT models. *Clin Pharmacol Ther* 100:170–178. doi:10.1002/cpt.361
- Blinova K, et al. (2018) International multisite study of human-induced pluripotent stem cell-derived cardiomyocytes for drug proarrhythmic potential assessment. *Cell Reports* 24:3582–3592. doi:10.1016/j.celrep.2018.08.079
- Blinova K, et al. (2019) Clinical trial in a dish: personalized stem cell-derived cardiomyocyte assay compared with clinical trial results for two QT-prolonging drugs. *Clin Transl Sci* 12:687–697. doi:10.1111/cts.12674
- Colatsky T, et al. (2016) The Comprehensive in Vitro Proarrhythmia Assay (CiPA) initiative: update on progress. *J Pharmacol Toxicol Methods* 81:15–20. doi:10.1016/j.vascn.2016.06.002
- Crumb WJ, et al. (2016) An evaluation of 30 clinical drugs against the comprehensive in vitro proarrhythmia assay (CiPA) proposed ion channel panel. *J Pharmacol Toxicol Methods* 81:251–262. doi:10.1016/j.vascn.2016.03.009
- Dutta S, et al. (2017) Optimization of an in silico cardiac cell model for proarrhythmia risk assessment. *Front Physiol* 8:616. doi:10.3389/fphys.2017.00616
- Fermini B, et al. (2016) A new perspective in the field of cardiac safety testing through the comprehensive in vitro proarrhythmia assay paradigm. *SLAS Discov* 21:1–11. doi:10.1177/1087057115594589
- Gintant G, et al. (2022) Assessing the fidelity of translation of non-clinical assays: a Pharma perspective. *Br J Pharmacol*. doi:10.1111/bph.15796 (read via abstract only)
- Guerrelli D, et al. (2024) hiPSC-CM electrophysiology: impact of temporal changes and study parameters on experimental reproducibility. *Am J Physiol Heart Circ Physiol* 327:H12–H27. doi:10.1152/ajpheart.00631.2023
- Hnatkova K, et al. (2019) Heart rate correction of the J-to-Tpeak interval. *Sci Rep* 9. doi:10.1038/s41598-019-51491-4
- Johannesen L, et al. (2016) Automated algorithm for J-Tpeak and Tpeak-Tend assessment of drug-induced proarrhythmia risk. *PLoS ONE* 11:e0166925. doi:10.1371/journal.pone.0166925
- Li Z, et al. (2017) Improving the in silico assessment of proarrhythmia risk by combining hERG channel-drug binding kinetics and multichannel pharmacology. *Circ Arrhythm Electrophysiol* 10. doi:10.1161/CIRCEP.116.004628
- Li Z, et al. (2019) Assessment of an in silico mechanistic model for proarrhythmia risk prediction under the CiPA initiative. *Clin Pharmacol Ther* 105:466–475. doi:10.1002/cpt.1184
- Malik M, et al. (2019) Implications of individual QT/RR profiles, part 1. *Drug Saf* 42:401–414. doi:10.1007/s40264-018-0736-1
- Park E, et al. (2018) Can non-clinical repolarization assays predict the results of clinical thorough QT studies? Results from a research consortium. *Br J Pharmacol*. doi:10.1111/bph.14101 (read via abstract only)
- Strauss DG, et al. (2021) Translational models and tools to reduce clinical trials and improve regulatory decision making for QTc and proarrhythmia risk (ICH E14/S7B updates). *Clin Pharmacol Ther*. doi:10.1002/cpt.2137 (read via abstract only)
- Vicente J, et al. (2018) Mechanistic model-informed proarrhythmic risk assessment of drugs. *Clin Pharmacol Ther* 103:54–66. doi:10.1002/cpt.896
- Goldberger AL, et al. (2000) PhysioBank, PhysioToolkit, and PhysioNet. *Circulation* 101:e215. doi:10.1161/01.cir.101.23.e215
- Johannesen L, et al. (2014) Differentiating drug-induced multichannel block on the electrocardiogram: randomized study of dofetilide, quinidine, ranolazine, and verapamil. *Clin Pharmacol Ther* 96:549–558. doi:10.1038/clpt.2014.155
- Johannesen L, et al. (2015) Late sodium current block for drug-induced long QT syndrome: results from a prospective clinical trial. *Clin Pharmacol Ther* 99:214–223. doi:10.1002/cpt.205
- Millard D, et al. (2018) Cross-site reliability of human induced pluripotent stem cell-derived cardiomyocyte based safety assays using microelectrode arrays: results from a blinded CiPA pilot study. *Toxicol Sci* 164:550–562. doi:10.1093/toxsci/kfy110
- O'Hara T, et al. (2011) Simulation of the undiseased human cardiac ventricular action potential: model formulation and experimental validation. *PLoS Comput Biol* 7:e1002061. doi:10.1371/journal.pcbi.1002061
- Pritchard JF, et al. (1985) Plasma protein binding of bepridil. *J Clin Pharmacol* 25:347–353. doi:10.1002/j.1552-4604.1985.tb02853.x
- Redfern WS, et al. (2003) Relationships between preclinical cardiac electrophysiology, clinical QT interval prolongation and torsade de pointes for a broad range of drugs: evidence for a provisional safety margin in drug development. *Cardiovasc Res* 58:32–45. doi:10.1016/S0008-6363(02)00846-5 (closed access; one excerpt read via a search result only)
- Sager PT, et al. (2014) Rechanneling the cardiac proarrhythmia safety paradigm: a meeting report from the Cardiac Safety Research Consortium. *Am Heart J* 167:292–300. doi:10.1016/j.ahj.2013.11.004
- Schocken D, et al. (2018) Comparative analysis of media effects on human induced pluripotent stem cell-derived cardiomyocytes in proarrhythmia risk assessment. *J Pharmacol Toxicol Methods* 90:39–47. doi:10.1016/j.vascn.2017.11.002
- Shinozawa T, et al. (2017) Recapitulation of clinical individual susceptibility to drug-induced QT prolongation in healthy subjects using iPSC-derived cardiomyocytes. *Stem Cell Reports* 8:226–234. doi:10.1016/j.stemcr.2016.12.014
- Smith AH, et al. (2007) Autonomic tone attenuates drug-induced QT prolongation. *J Cardiovasc Electrophysiol* 18:960–964. doi:10.1111/j.1540-8167.2007.00901.x
- Stillitano F, et al. (2017) Modeling susceptibility to drug-induced long QT with a panel of subject-specific induced pluripotent stem cells. *eLife* 6. doi:10.7554/elife.19406
- Strodthoff N, et al. (2023) PTB-XL+, a comprehensive electrocardiographic feature dataset. *Sci Data* 10. doi:10.1038/s41597-023-02153-8
- Vicente J, et al. (2015) Comprehensive T wave morphology assessment in a randomized clinical study of dofetilide, quinidine, ranolazine, and verapamil. *J Am Heart Assoc* 4. doi:10.1161/JAHA.114.001615
- Vicente J, et al. (2016) Electrocardiographic biomarkers for detection of drug-induced late sodium current block. *PLoS ONE* 11:e0163619. doi:10.1371/journal.pone.0163619
- Vicente J, et al. (2019) Assessment of multi-ion channel block in a phase I randomized study design: results of the CiPA phase I ECG biomarker validation study. *Clin Pharmacol Ther* 105:943–953. doi:10.1002/cpt.1303
