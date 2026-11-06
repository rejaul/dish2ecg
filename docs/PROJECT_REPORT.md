# From the dish to the ECG

*A public-data comparison of stem-cell heart cells and clinical QT results*

Rejaul Karim, Campbellsville University, Campbellsville, Kentucky, USA

dish2ecg version 0.2.2 · 6 November 2026 · https://github.com/rejaul/dish2ecg

## Summary

Drug safety teams now screen new compounds for torsades risk partly in a dish, using heart cells grown from human induced pluripotent stem cells (hiPSC-CMs). The dish gives a field potential duration (FPD). It looks like the QT interval on a patient's ECG, but it is a different measurement. Published work often puts the two side by side without saying how each was corrected or on what concentration scale. We wanted to know how far that comparison can be tested with public data alone.

To do that we wrote an open Python package, dish2ecg. It sends the dish data and the clinical ECG data through the same correction steps and compares them as multiples of the free drug concentration that patients reach. Before relying on each stage we checked it against a published analysis. The clinical side matched published tables in 8 of 8, 5 of 6 and 6 of 8 comparisons, and our dish-side heart-rate correction fell within 5 ms of an independent implementation in 17 of 18 studies.

Only four drugs have both dish data and clinical concentration data: dofetilide, quinidine, ranolazine and verapamil. At the same free concentration the dish gets the direction right for three of them and overstates the effect by a factor of 1.8 to 3.5. For verapamil it goes the wrong way: +3.5 ms in volunteers and −86 ms in the dish. Sorting 28 drugs by torsades risk is an easier job. On the 16 drugs that the FDA's in-silico score (qNet) was not fitted on, the dish ranks drugs about as well as qNet does (Spearman correlation with risk class 0.81 against 0.84) and is as good at clearing low-risk drugs. qNet is better at flagging the high-risk ones.

Four problems remain open. Public data cannot separate the type of dish measurement from the laboratory and the culture medium. The clinical recordings are only 10 seconds long, so the usual correction for QT "memory" cannot be applied. Two clinical studies disagree by 6 to 7 ms for the same drugs and we cannot explain it. And four drugs are too few to calibrate anything. We also withdrew four of our own claims after testing them. They are listed in section 5.

Readers who want the physiology first can start at section 12. The reading list is in section 17.

## 1. The question

A cardiac safety package has to answer one question twice. Does the compound delay repolarisation of the ventricle, and by how much, at the concentration a patient will have in the blood? In the clinic the answer comes from the QT interval of a 12-lead ECG, corrected for heart rate (QTc). In the laboratory, since the CiPA initiative (Sager 2014; Colatsky 2016; Fermini 2016), part of the answer comes from a plate of hiPSC-CMs sitting on electrodes. The plate records a field potential, and its duration (FPD) is corrected for beat rate in much the same way.

The two readouts are often placed side by side. A reader can rarely check how each was corrected, which concentration scale was used, or whether the dish number is being read as if it were milliseconds of QT. The largest public dish dataset is the CiPA myocyte study of Blinova and colleagues (2018): 28 drugs tested at 10 sites on two commercial cell lines. The largest public source of clinical concentration-QT data is PhysioNet (Goldberger 2000), which holds three early-phase studies from the FDA and its collaborators.

Our question was whether both can be put through one correction algebra and compared on one concentration scale, and how far public data let us test that. We did not set out to show that the dish predicts QT. We made one rule at the start and the code enforces it: a dish endpoint is never expressed as a QT interval in milliseconds. The package raises an error if asked to.

*Figure 1. What the two dish readouts measure. (a) The field potential from an electrode array. FPD runs from the depolarisation spike to the peak of the repolarisation wave, which mainly reflects the potassium current IKr. (b) The action potential from a voltage-sensitive optical probe. APD90 is the time to 90 % repolarisation. Neither is a QT interval. (This figure is in the HTML and PDF editions of the report.)*

## 2. Data

Everything came from public sources. We did not redistribute any third-party data. The repository holds code, derived results and pointers to the originals.

| Source | What it holds | Size |
|---|---|---|
| CiPA myocyte study, Blinova et al. 2018 (workbook supplied by the project owner) | Well-level corrected change in field potential duration (ddFPDc), arrhythmia flag, 4 concentrations per drug | 28 drugs, 8,859 wells, 10 sites, 6 platforms, 2 cell lines (iCell², Cor.4U), 15 dataset-by-platform combinations |
| ECGRDVQ (PhysioNet) | Dofetilide, quinidine, ranolazine, verapamil and placebo: intervals, plasma levels, waveforms | 22 subjects, crossover |
| ECGDMMLD (PhysioNet) | Dofetilide with or without mexiletine or lidocaine; moxifloxacin with or without diltiazem | 22 subjects, crossover |
| ECGCIPA (PhysioNet) | Ranolazine, verapamil, lopinavir/ritonavir, chloroquine (parallel arms); dofetilide with or without diltiazem (crossover) | 60 subjects |
| Millard et al. 2018, Supplementary Table 1 | Beat period, FPD and corrected FPD for 18 studies | 18 rows |
| FDA CiPA repository (pinned commit 539ff80) | Ion-channel potency fits for 28 drugs; 2,000 bootstrap hERG fits per drug; the ORd-CiPA model | 28 drugs |
| Crumb et al. 2016 Table 1, product labels | Free Cmax; plasma protein binding | 30 drugs; 15 reference rows |

The clinical anchor is 11 molecules in 104 healthy volunteers. Nobody in these studies had an arrhythmia, which matters later. In the dish table, 2,067 of the 8,859 values are missing, and 1,887 of those are exactly the wells flagged as early afterdepolarisations (EADs), meaning wells that stopped beating normally.

## 3. Methods

### 3.1 One cascade for both sides

The idea is simple. A raw interval is corrected for rate. The subject's or well's own baseline is subtracted (Δ). The placebo or vehicle change is subtracted (ΔΔ). A concentration-response model is fitted, and the result is put on a shared concentration scale. Each step is recorded in a provenance tuple that travels with the number. The code refuses steps in the wrong order, refuses to subtract a control from something that has no baseline, and refuses to subtract twice.

![Figure 2. The processing pipeline, drawn. Rows are the shared cascade steps and columns are the three arms. Every panel below the top row shows this repository's own data: the 18 Millard studies moving from FPD to FPDc, 439 drug-free clinical records, the dofetilide baseline-and-placebo waterfall at 2.5 h, and the three fitted concentration-responses with 90 % bands. The ICH E14 10 ms line is drawn on the clinical panel, and the in-silico curve is cut off where the mid-myocardial cell stops repolarising. The top row is schematic, because no raw waveform is redistributed here. Red marks the two places the method refuses to go: expressing FPD as a QT interval, and treating a nominal dish concentration as a measured free one.](../results/pipeline_diagram.png)

### 3.2 Rate correction

We used the power form Ic = I / RR^α, with RR in seconds. Fridericia's α = 1/3 is the default. We fitted exponents on drug-free records only (baseline and placebo), so that the drug effect is not absorbed into the correction. They came out at 0.301 in ECGRDVQ and 0.299 in ECGDMMLD. That is what healthy volunteers should give, and it is a useful check on the fitting code. For J-Tpeak we used the published exponent 0.58 (Johannesen 2014) and confirmed it on 5,740 ECGCIPA records: a regression recovers 0.5799 with a residual SD of 0.15 ms. The exponent is not ours.

### 3.3 Concentration-response models

On the clinical side we fitted a linear mixed model of ΔΔQTcF against concentration, with a random intercept for each subject. The effect is read at the geometric-mean observed Cmax with a two-sided 90 % interval. Over 200 simulated crossover studies this interval covered the true effect 86.5 % of the time against a nominal 90 %, so our clinical intervals are a little too narrow.

On the dish side we used the same model with the dataset (site, cell line and platform) as the random effect, and log10 of the multiple of free Cmax as the predictor. We needed the log scale because tested concentrations span up to 3,500-fold, and a linear fit is dominated by its highest point. Wells that stop beating cannot be measured, and the loss rises steeply with dose: 0.9 %, 4.7 %, 26.8 % and 60.8 % across the four concentrations. So every dish fit is reported twice, once over all four concentrations and once over the lowest three. This bounds the problem. It does not fix it.

### 3.4 The bridge between dish and clinic

We converted plasma concentration from ng/mL to nM with the molecular weight, multiplied by the fraction unbound to get a free concentration, and divided by the free Cmax to get multiples of Cmax. Both sides were evaluated at the geometric-mean free Cmax that the volunteers actually reached, not at a nominal 1×. Quinidine reached only 0.29× and verapamil 0.50× of the reference value, so a nominal 1× would have meant extrapolating the clinical fit. Because the same reference scales both axes, the bridge does not depend on which published Cmax table is used.

One asymmetry survives all of this, and it is large. The clinical axis is a measured free plasma concentration. The dish axis is the nominal concentration added to the well. Section 7 comes back to it.

### 3.5 Waveforms

The PhysioNet studies include raw waveforms and median beats as well as interval tables. We decoded them (WFDB formats 16 and 32, with per-microvolt gains), detected R peaks, and re-measured the end of the T wave with a tangent rule. Within each 10-second strip we computed the spread of QTcF between the slowest and the fastest beat.

### 3.6 In-silico work

We ran the O'Hara-Rudy human ventricular model (O'Hara 2011) in myokit. The set was 408 simulations: the 12 CiPA training drugs, three cell types (endocardial, epicardial and mid-myocardial), and 11 or 12 concentrations per drug, with no failed runs. Four of the 12 drugs are the bridge drugs. For the risk comparison we used the FDA's ORd-CiPA v1.0 model (Dutta 2017) and computed qNet for all 28 drugs, each at four doses and for 20 bootstrap hERG parameter sets. That is 2,240 simulations, again with no failures. A coupled one-dimensional cable could not be run, for a software reason described in section 10.

### 3.7 Checking the work

The repository has 142 automatic tests, counted with pytest's collect-only mode. Every result file has a producing script or is labelled as made by hand. An audit near the end found 14 of 38 result files with no producer. Nine are now reproducible, one was deleted, and four are labelled ad hoc in `results/README.md`.

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

Two of these lines need a comment. The dofetilide line closed what had looked like our one unexplained gap. It turned out to be a disagreement between two published analyses of the same trial, and our value matches the primary paper. In the dish rate-correction check, study 9 is off by 22.0 ms and we do not know why. It stays in the test suite as an asserted discrepancy so it cannot be forgotten.

### 4.2 The four-drug bridge

![Figure 3. Change in repolarisation duration against multiples of free Cmax for the four drugs that exist in both data sets. Clinical ΔΔQTcF is blue and dish ΔΔFPDc is orange, with the matched-exposure readout in each panel. The dish agrees in sign for three drugs and reverses for verapamil.](../results/dish_vs_clinic.png)

| Drug | Free Cmax multiple | Clinic ΔΔQTcF (ms) | Dish ΔΔFPDc (ms) | Dish ÷ clinic |
|---|---|---|---|---|
| Dofetilide | 1.07 | +73.1 [68.5, 77.7] | +132.0 [93.1, 170.8] | 1.81 |
| Quinidine | 0.29 | +74.9 [69.3, 80.5] | +135.1 [82.2, 187.9] | 1.80 |
| Ranolazine | 0.96 | +11.2 [7.6, 14.8] | +38.8 [22.8, 54.7] | 3.46 |
| Verapamil | 0.50 | +3.5 [−0.1, 7.0] | −86.3 [−98.9, −73.6] | −24.7 |

The two pure hERG blockers are inflated by almost the same factor (1.81 and 1.80), which looks like a property of the assay and not of the drugs. The dish intervals are several times wider than the clinical ones, even with 128 to 278 wells per drug. The reason is variance between sites, platforms and cell lines, not noise between wells. Four drugs are the whole overlap. A search of PhysioNet, Zenodo and the Scientific Data descriptors found no other public well-level field-potential dataset, and no other public single-agent concentration-QTc set.

### 4.3 Why verapamil reverses

At the Crumb free Cmax, dofetilide blocks hERG by 24 % and has little effect on the calcium and sodium currents. Verapamil blocks hERG by 15 % and the L-type calcium current by 16 %. The usual reading is that in volunteers the calcium block offsets the potassium block, so QTc barely moves. In the dish, the likely reason for the reversal is that immature hiPSC-CMs have little of the inward rectifier IK1 and beat spontaneously. Calcium current then carries more of their plateau, and blocking it shortens the field potential. We did not test this.

We did test whether verapamil's own channel profile is enough to cause the reversal, by running it through the ORd model. It is not. The model prolongs the action potential slightly for verapamil, as the clinic does, and it gets the sign right for all four drugs. So the dish is the part that is off. The model is far from perfect on size. Ranolazine comes out at about +37 ms in the mid-myocardial cell against +11 ms in volunteers, which is close to the dish's +39 ms. We read the in-silico arm as a check on direction and not as a prediction of QTc.

The uncoupled model also says something about the QT surrogate. The mid-myocardial (M) cell had the longest APD90 in 121 of 136 conditions. In all 15 exceptions, which sit at 3.78× free Cmax or above, the M cell itself had failed to repolarise (section 10). So in every condition where it repolarised normally it was the longest, and an uncoupled "QT" is simply the longest cell, with an implied dispersion of about 102 ms where a real ventricular wall has 30 to 40 ms. The common mapping that takes a weighted average of action potential durations has no mechanism behind it. The control APD90 values were 265.5 ms (endocardial), 228.0 ms (epicardial) and 329.9 ms (mid-myocardial).

### 4.4 Endpoint, laboratory and medium cannot be separated

Blinova 2018 stores two different measurements in one column, `ddFPDc`. Thirteen of the 15 datasets are electrode FPD. The other two are the same optical platform on the two cell lines (15.2 % of wells), where the measurement is APD90 from a voltage-sensitive probe. The paper calls the two endpoints "equivalent measures". As a scale they are not. The optical arm reads 1.44 times the consensus of the other nine sites (interval 1.15 to 1.72), the highest of the ten. In the paired comparison, optical = 1.507 × electrode + 24.8 ms (n = 97).

The difficulty is that the optical platform was used at one site only, and that site is also the only one that used serum-free medium. Schocken and colleagues (2018) showed that medium alone moved baseline corrected FPD by 104 ms in iPSC-CMs (575 against 471 ms). So endpoint, site and medium move together, and no public dataset separates them. We can say that the optical arm reads higher. We cannot say why.

*Figure 4. Three rungs from cell to patient. Left: gain of each site relative to the other nine, with the single optical site highest. Middle: optical against electrode readings of the same drug and concentration, slope 1.51. Right: clinical, dish and in-silico effects for the four bridge drugs. (This figure is in the HTML and PDF editions of the report.)*

### 4.5 The clinical side has its own floor

Each clinical ECG record in the waveform data is 10 seconds long. QT takes one to two minutes to settle after a change in heart rate, and 10 seconds is only 8 to 17 % of one adaptation time constant. Most of the RR history that sets the QT in a strip happened before the strip began. This makes the hysteresis correction unidentifiable here, which is a stronger statement than "not implemented".

The cost is measurable. Across 238 records, QTcF computed on the slowest and the fastest beat within one strip differs by a median of 12.6 ms, and by more than 10 ms in 66.8 % of records. The studies take replicate ECGs at each time point, which reduces the effect on a mean. Even so, it shows how little a single 10-second strip fixes about steady-state QT.

The end of the T wave is also rate dependent. Our tangent rule differs from the study's own markings by a median of −4.5 ms in ECGRDVQ and −7.1 ms in ECGCIPA, and the difference changes with RR (−0.032 ms per ms of RR, p = 0.026). A constant offset between two delineation methods cancels in a double-delta analysis. A rate-dependent offset does not, if heart rate differs between baseline and treatment. PTB-XL+, where two commercial algorithms measured the same 21,795 ECGs, shows the same sign, with a difference of about −5.9 ms across a 300 ms change in RR.

![Figure 5. Waveform findings from 238 clinical records. Left: beat-to-beat RR inside single ten-second recordings. Middle: spread of QTcF between each window's slowest and fastest beat. Right: difference between our tangent T-end and the study's, against RR.](../results/waveform_findings.png)

### 4.6 Sorting 28 drugs by risk

The assay exists to sort drugs into torsades risk classes (8 high, 11 intermediate and 9 low in the CiPA coding). We compared dish endpoints against the in-silico qNet score on all 28 drugs. Twelve of them were used to train qNet and 16 were held out.

We first tried this on the 12 training drugs only and concluded that the dish "detects danger but does not grade safety". That conclusion was wrong. It came from looking at the drugs on which qNet had been tuned, with only 12 points. Once the FDA's bootstrap inputs for all 28 drugs turned out to be public, we could test the 16 held-out drugs. Our qNet reproduces the published point fit on 11 of the 12 training drugs to within 0.0006 µC/µF.

| On the 16 held-out drugs | Spearman ρ with risk class | AUC, high vs rest | AUC, low vs rest |
|---|---|---|---|
| qNet (in silico) | 0.840 [0.52, 0.95] | 1.000 | 0.909 |
| ddFPDc at Cmax (dish) | 0.809 [0.55, 0.93] | 0.917 | 0.945 |
| max ddFPDc (dish) | 0.722 [0.38, 0.90] | 0.812 | 0.945 |
| Arrhythmia-event fraction (dish) | 0.717 [0.37, 0.90] | 0.875 | 0.891 |

qNet separates every high-risk drug from every other drug, on all 28. That is the only comparison in the table with no overlap. On ranking and on clearing low-risk drugs, the dish and qNet cannot be told apart at n = 16, because the intervals overlap almost completely. qNet's two false alarms, both held out, are metoprolol (0.0702) and tamoxifen (0.0781).

*Figure 6. (a) qNet for 28 drugs, median over 20 bootstrap hERG fits. Open circles are the 16 held-out drugs, and every high-risk drug falls to the left of the dashed line. (b) Areas under the ROC curve on the 16 held-out drugs. (This figure is in the HTML and PDF editions of the report.)*

Two qNet traps cost us time and may cost someone else. The current set that qNet integrates excludes the sodium current, and the L-type calcium current acts as the calcium carrier only. Including INa gave 0.0044 where the correct value is 0.0785. The trapping constant Kt = 3.5e-5 is fixed in CiPA but defaults to zero in the model file, which switches trapping off. With Kt at zero, dofetilide at 4× prolonged APD90 by 28 ms instead of 130 ms.

### 4.7 A 6 to 7 ms shift we cannot explain

Ranolazine, verapamil and dofetilide were measured in both ECGRDVQ and ECGCIPA. At a matched plasma concentration ECGCIPA reads higher every time: by 8.2 ms for ranolazine, 5.0 for verapamil and 9.7 for dofetilide. ECGCIPA subtracts one shared placebo group from every arm. That accounts for roughly 2 to 4 ms and for essentially all of the verapamil shift. For ranolazine and dofetilide, 6 to 7 ms remains. With no placebo subtraction at all the shifts are 6.0 and 7.0 ms.

We tested the most obvious explanation, the T-end rule, on 3,163 median beats from ECGRDVQ and 2,515 from ECGCIPA. It moves the drug effect by −0.94 and −1.83 ms in ECGRDVQ (ranolazine and dofetilide) and by −0.61 and +0.25 ms in ECGCIPA. The difference between the studies is −0.32 ms for ranolazine and −2.08 ms for dofetilide. That is at most about a third of the dofetilide shift and a twentieth of the ranolazine shift. The placebo-corrected change in RR is only −33 to +18 ms, so a rate-dependent rule would predict a movement of just +0.3 to +1.0 ms. The shift stays open. We had briefly written the T-end rule up as a live contributor, and that claim is withdrawn.

### 4.8 Free fraction for the other CiPA drugs

The bridge needs the fraction unbound (fu) for each drug. For the four bridge drugs it is fixed from the literature: dofetilide 0.35, quinidine 0.16, ranolazine 0.38, verapamil 0.10. Four other CiPA drugs had no fu in our table, so we looked them up late in the project.

| Drug | fu (point) | Range | Basis |
|---|---|---|---|
| Bepridil | 0.005 | 0.003–0.010 | Label: more than 99 % bound. Equilibrium dialysis (Pritchard 1985): more than 99.7 % bound. Binding saturates. |
| Chlorpromazine | 0.055 | 0.03–0.08 | 92–97 % bound (drug monograph, not an FDA label) |
| Cisapride | 0.0225 | 0.020–0.025 | Label: 97.5–98 % bound |
| Terfenadine | 0.03 | none given | Label (1995): 97 % bound |

None of the four has a current FDA label, so the text comes from reproductions on drug-information sites. Bepridil's point value is a bounded estimate and not a measurement, because its binding depends on concentration. These four drugs have no patient data and are not in the bridge, so no result above depends on them.

The search did turn up one thing worth keeping. Terfenadine's free Cmax is 4 nM in the CiPA table and 0.3 nM in Crumb and Blinova, a 13-fold difference. With fu = 0.03 these imply total concentrations of about 63 and 4.7 ng/mL, against a reported normal-dose peak near 1.5 ng/mL. An excerpt of Redfern 2003, which we saw only as a search snippet and not in full, says the 4 nM value is a concentration reached with a P450 inhibitor present. If so, the two tables describe different exposures, and a "multiples of Cmax" axis means something different for terfenadine depending on which one is used. Crumb 2016 and Li 2017 do not tabulate protein binding, contrary to what an earlier note in this repository said.

## 5. What we got wrong

Four claims were made and then withdrawn after a test went against them.

1. We said that correcting for free drug concentration would widen the gap between dish and clinic. Schocken 2018 found that serum-containing medium made 8 of 18 drugs arrhythmic at lower concentrations, the opposite of what protein binding alone predicts. The sign is not established.
2. We blamed unpaired placebo subtraction in a parallel design for the difference between two clinical studies. That did not survive a direct test.
3. We said the dish detects danger but does not grade safety. This was an artefact of using only the 12 training drugs.
4. We said the T-end rule was a live contributor to the 6 to 7 ms shift. The beat-level test showed it is not.

We also fixed smaller errors in the repository's documents. A median of 0.6 ms for the Millard comparison was really 0.65 ms. A paper had been credited with naming one "best" predictor when it names three. A test count no longer matched the suite. An early note claimed that two papers tabulate plasma binding when they do not. The confidence intervals in the held-out risk table had gone stale, so that the qNet Spearman interval read [0.55, 0.95] instead of [0.52, 0.95]. And a 2019 heart-rate-correction paper was attributed to the wrong first author in the citation file.

## 6. Discussion

The dish and the clinic agree on direction for most drugs and disagree on size. That is not surprising. A denervated, immature, single-layer preparation with no His-Purkinje system should over-read a hERG blocker. A human study found that autonomic tone attenuates drug-induced QT prolongation. Ibutilide prolonged QTc by 43 ± 8 ms with tone intact, and by a further 68 ± 5 ms after double autonomic blockade, on top of the 16.7 ± 9.4 ms that the blockade caused by itself (Smith 2007). That fits the direction of our 1.8 to 3.5 fold inflation. It is a direction and not an explanation, and we claim no mechanism for the size.

The more useful result may be that the dish does well on a task closer to what it is used for. A drug-level sort into risk classes does not need FPD in milliseconds, and on that task the dish at its Cmax is about as good as the FDA's model at identifying low-risk drugs. The four-drug bridge asks the dish for a quantity it was never designed to give. We do not think the two results conflict. They suggest the assay is better as a ranking tool than as a measuring tool.

On individual prediction the evidence is weaker still. In the one study that grew cells from the same volunteers who took the drugs, subject-level dish slopes did not correlate with subject-level clinical slopes (Spearman ρ = −0.11, p = 0.70 for dofetilide). The same subjects' clinical responses to the two drugs did not correlate with each other either (ρ = 0.10, p = 0.73), so part of the weakness lies on the clinical side. Earlier studies reported positive correlations in other settings (Stillitano 2017; Shinozawa 2017). We use the bridge at the drug level only.

What is new here is modest and mostly about method. We found no open implementation that runs the same correction cascade on dish and clinic data, carries provenance, and refuses to turn FPD into QT. We also found no public side-by-side at free-Cmax multiples that re-derives the clinical side from open ECG data, though effect-at-Cmax in the dish is not new. The measurement of informative missingness (1,887 of 2,067 nulls are EAD wells) is a quantification and not a discovery. The J-Tpeak exponent is not novel. These verdicts rest on a limited literature screen, and earlier screens in this project were overturned by sources the first pass missed, so none of them is a claim of priority.

## 7. Limits

The bridge rests on four drugs, so every statement about its size carries wide uncertainty and none of it supports a calibration. The clinical population is healthy adults aged 19 to 35 with no arrhythmic events. That is the group least likely to show the outcome the dish is supposed to anticipate.

The dish concentration axis is nominal, and this is the weakest link in the comparison. What is pipetted into a well is not what the ion channels see. Drug sticks to the plate plastic and to any tubing, and lipophilic compounds partition into the cell membrane and into intracellular lipid. This loss is separate from binding to serum protein in the medium. It depends on the compound, and because binding sites are finite the fractional loss is largest at the lowest concentrations, which is where a careful assay works. The patch-clamp literature calls it nonspecific binding or drug loss (Alvarez Baron 2025).

The size can be large enough to change a conclusion, and it varies widely by compound. A study that verified concentrations by LC-MS/MS in the recording rig, instead of trusting the nominal value, found hERG-experiment losses of 6 % for ranolazine, 19 % for chloroquine, 41 % for verapamil and 58 % for vanoxerine. Accounting for the loss moved the hERG IC50 to 2.2 nM for vanoxerine and to 0.75 µM for chloroquine. The same paper notes that earlier studies all plotted nominal concentrations because none had measured the in-vitro concentration at all.

Two of those four compounds, ranolazine and verapamil, are bridge drugs here, and their losses differ sevenfold. A per-compound loss does not cancel when two drugs are compared on a shared exposure axis. We cannot carry the numbers across, because those were patch-clamp rigs with flowing perfusion and our data come from static MEA wells, so the mix of mechanisms differs: tubing and perfusion lines there, plate plastic and cell lipid here. What carries over is that the loss is real, depends on the compound, and is large enough to move a potency estimate.

If free concentration in a well is below nominal, every dish curve belongs further to the left on the exposure axis. The cells would then be more potent per free molecule than plotted, and the dish's over-read at matched exposure would be larger. That is the direction of one mechanism taken alone, not a net correction. A related claim about protein binding was made and withdrawn earlier in this project (section 5, item 1), because the one experiment that bears on the net sign ran the other way: 8 of 18 drugs became arrhythmic at lower nominal concentrations when serum was present (Schocken 2018). Nobody measured free concentration in these wells, so the size and even the net sign of the full correction are open. We are not reinstating the withdrawn claim.

The in-silico comparison is out of sample for the 16 held-out drugs, but qNet was developed with the 12 training drugs in view, and with n = 16 the intervals are wide. The coupled one-dimensional cable did not run (section 10).

Some technical limits are smaller. The clinical model under-covers (0.865 against 0.90). We did not compare our degrees-of-freedom rule with Kenward-Roger. A linear model fitted to verapamil's curved response caused both misses against Vicente 2019. Residual diagnostics were not run.

## 8. What would settle the open questions

1. Endpoint, site and medium. A study in which the same cells and drugs are read by both electrodes and optical probes, in the same medium, at the same site.
2. Hysteresis. Beat-level or longer ECG recordings with the RR history, which the public tables do not carry.
3. The 6 to 7 ms shift. Refit both studies with the full white-paper model (treatment indicator, nominal-time effects, baseline covariate), refit ECGCIPA day by day, and compare extraction methods, ideally on the waveforms of both.
4. Missing wells. A tipping-point analysis that gives each missing well a duration equal to the observed maximum plus δ, and reports the δ at which a conclusion changes. This is not done.
5. A usable bridge. More single-agent clinical concentration-QTc data that are open. None exist publicly beyond these four.
6. The cable model. A computer with an OpenCL driver that exposes double precision at the platform level.
7. The dish exposure axis. Measurement of free compound in the well at the time of recording, by LC-MS/MS on the medium, for each drug and each site. Until that exists the dish axis stays nominal, and the size of the dish-clinic gap cannot be read as a property of the cells alone.

## 9. What this work is, and is not

I did not grow cells, record a field potential or run a clinical study. All the dish data come from the CiPA myocyte study and all the clinical data from PhysioNet. The contribution is the curation and modelling half. I assembled heterogeneous public data into an auditable form, implemented the correction cascade with provenance, fitted mixed-effects models that carry site and protocol variability explicitly, validated each stage against independent published analyses, decoded the clinical waveforms, ran the in-silico arm, and tested my own claims. Where the data could not answer a question, the result is a precise statement of why.

## 10. Challenges we met

These are practical obstacles, as opposed to the scientific limits in section 7. Each one shaped a result.

### 10.1 The dish data arrive already corrected

The published workbook ships ΔΔFPDc, not raw field potential durations and beat periods. Baseline, vehicle and rate correction were applied by the study, so our cascade cannot be run on the dish data. We verified the numbers instead (r = 0.999993 against all 420 published panels) and ran the same code on an independent dataset (Millard 2018, 17 of 18 within 5 ms). The consequence is that the in-vitro baseline and vehicle steps are tested only on synthetic data and on the clinical corpus, never on real raw wells.

### 10.2 Two endpoints in one column

`ddFPDc` holds electrode FPD for 13 datasets and optical APD90 for 2. The paper calls them equivalent measures. As a scale they are not: the optical arm reads 1.44 times the consensus of the other nine sites. But the optical platform was used at one site, which is also the only serum-free site, and medium alone shifts baseline FPDc by 104 ms (Schocken 2018). Endpoint, site and medium move together, and no public dataset separates them.

### 10.3 Missing wells

Of the 8,859 values, 2,067 are missing, and 1,887 of those are exactly the wells flagged as arrhythmic. The loss rises from 0.9 % to 60.8 % across the four concentrations. The durations of those wells do not exist in the published data, so we bound the problem by fitting everything twice and do not impute. A tipping-point analysis is the obvious next step. It has not been done.

### 10.4 Small defects in the source tables

Fifteen wells at the top concentration are under-coded for diltiazem, nitrendipine and verapamil, all at one site on one platform. One cell (clozapine, site 1, iCell², top concentration) has five measurable wells but no plotted value. Neither changes a conclusion. Both are recorded in the validation report and not quietly dropped.

### 10.5 A validation that still fails

Study 9 (iCell² on CardioECR) is 22.0 ms off in the rate-correction check, with a beat-period coefficient of variation of only 0.061, so averaging order does not explain it. We could not account for it. It is asserted in a test so it cannot be forgotten.

### 10.6 A quality flag that misses the failure

In the in-silico arm the mid-myocardial cell stops repolarising at high exposure, from about 3.8 times free Cmax upward. Its APD90 is then either missing or measured from a depolarised resting potential, which gives values like −272 ms. The run's own `ok` flag is true for all 408 runs, so it cannot be trusted, and the physiology has to be checked directly (a resting potential below −60 mV and a finite APD90). We found this while drawing Figure 2. None of the affected rows was used in a reported result.

### 10.7 Technical blockers

The coupled one-dimensional cable did not run. The OpenCL driver on our remote machine exposes double precision on the device but not on the platform, so the solver fell back to single precision and the model returned NaN. A single-cell model file also lacks the diffusion-current binding that the conductance setter needs, and in that state the call silently does nothing. Two papers behind paywalls could not be read, and one number in this report rests on an excerpt rather than a full text (Redfern 2003).

### 10.8 Published sources that disagree

One published table disagrees with its own study's primary paper on the dofetilide effect. Our value matches the primary paper to 0.2 ms. The two reference tables for free Cmax disagree for terfenadine by a factor of 13, which we traced to one of them being an exposure measured with a P450 inhibitor present.

## 11. Questions for the study authors

Several of the limits above come less from the public data than from what the data release records. The group that ran the CiPA myocyte study could answer most of them directly. The questions are ordered by how much they would change what can be concluded.

### The pooled endpoint

Q1. Did any site run both the electrode and the optical platform on the same cells, in the same medium? Even a small paired set would separate endpoint from site and medium, which no public dataset does.

Q2. Was the 1.44-fold higher reading at the optical site noticed during the study? Is it understood as an endpoint effect, a medium effect or a site effect?

Q3. On what basis were FPD and APD90 treated as equivalent measures: agreement in rank, in drug effect, or in absolute scale? Our paired comparison gives optical = 1.51 × electrode + 24.8 ms.

### The arrested wells

Q4. For a well that develops an arrhythmic event, is any duration recorded before the event, or is the measurement abandoned at that point?

Q5. Is there a well-level record of why each value is missing, separating arrhythmia from signal loss, from a well that never beat, from exclusion by a quality rule?

Q6. How were the missing wells handled in the published concentration-response analysis?

### Exposure

Q7. Was the concentration in the well ever verified analytically, for any drug or any site? Recent patch-clamp work found losses of 6 % to 58 % between nominal and measured concentration.

Q8. Was anything done to limit adsorption to the plate and the tubing, and did the protocol differ between the serum-free site and the others?

Q9. Why does the free Cmax used in this study differ from the CiPA table for terfenadine by a factor of 13, and was the choice deliberate?

### Protocol

Q10. How long after drug addition was each reading taken, and did that interval differ by site? One published study reports beat rate rising 22.6 % and FPD falling 7.7 % over a 20-minute equilibration, which is the same order as the drug effects being measured.

Q11. What exactly do the baseline and vehicle corrections subtract: the same well before dosing, the plate mean, or the site mean, and over what time window?

Q12. Are the raw, uncorrected FPD and beat-period values kept anywhere? They would let the correction cascade be validated on real wells and not only on synthetic data.

### Specific results

Q13. Does the shortening of FPD by verapamil appear at every site, and was it expected from its calcium-channel block? It is the one drug where the dish and the clinic disagree in sign, and a human ventricular model given verapamil's own channel profile agrees with the clinic.

Q14. Is anything known about study 9 (iCell² on CardioECR) that would explain a 22 ms discrepancy in a rate-correction check that the other 17 studies pass?

Q15. After the subject-matched study of 2019, which found no within-person correlation, is a larger one planned? Our reading is that the clinical side contributed to that null, since the same subjects' responses to two drugs did not correlate with each other either.

## 12. The biology underneath

This section is background. A reader who wants the physiology before the methods can start here.

### 12.1 The action potential, and why the dish's is different

A ventricular cell at rest sits near −85 mV, held there by the inward rectifier potassium current IK1. A stimulus opens fast sodium channels. Sodium floods in and the cell depolarises within a millisecond or two (phase 0). Calcium then enters through L-type channels (ICaL) and holds the cell depolarised in a plateau, while potassium currents push the other way. The plateau ends when the potassium currents win. These are the rapid delayed rectifier IKr, carried by the hERG channel, and the slow one, IKs. The cell returns to rest (phase 3) and IK1 holds it there (phase 4). APD90 is the time from the upstroke until the cell has recovered 90 % of the way back.

hiPSC-CMs are immature. They have little IK1, so they rest depolarised near −72 mV rather than −85 mV, and they beat spontaneously instead of waiting to be stimulated. That is not only a defect. It is a prerequisite for the assay, because a quiescent adult myocyte would have to be paced. But it changes the balance of currents, and that is the root of most of what follows.

### 12.2 What each blocker does

The table uses the measured block at each drug's free Cmax (`results/channel_block_at_crumb_free_cmax.csv`).

| Drug | hERG (IKr) | ICaL | Late sodium | Net effect |
|---|---|---|---|---|
| Dofetilide | 24.0 % | 0.4 % | 3.5 % | Almost pure IKr block, so prolongation is unopposed |
| Quinidine | 46.7 % | 8.1 % | 3.8 % | IKr block dominant, with Ito 13.9 % and IKs 8.3 % |
| Ranolazine | 22.6 % | none | 21.1 % | IKr block partly offset by late-sodium block |
| Verapamil | 14.6 % | 16.2 % | 0.5 % | IKr block offset by a comparable calcium block |

Blocking IKr removes repolarising current, so the plateau lasts longer and the cell takes longer to recover. Blocking ICaL removes depolarising current, which shortens the plateau. Blocking late sodium current does the same by a different route. A drug that does both, like verapamil, ranolazine or mexiletine, has two effects that partly cancel, and which one wins depends on the cell.

This is why verapamil is such a useful failure. In an adult ventricular myocyte the two effects roughly balance and QTc barely moves (+3.5 ms in volunteers). In an immature cell with little IK1, calcium current carries proportionally more of the plateau, so removing it shortens the field potential sharply (−86 ms in the dish). When verapamil's own channel profile is run through a human ventricular model, the model prolongs slightly, like the clinic. The inversion is therefore a property of the preparation and not of the compound.

### 12.3 From one cell to a field potential

A multielectrode array does not record the transmembrane potential. The electrodes sit outside the cells, under a sheet of them, and record the extracellular field that the sheet generates. That field is closer to the time derivative of the membrane potential than to the potential itself. The sharp downward spike is the synchronised sodium upstroke. The slow hump much later is the repolarisation wave, and its peak is the conventional end of the field potential duration. FPD therefore follows repolarisation and is sensitive to IKr block, but it is not an action potential duration, and the two are not interchangeable (section 4.4 and `results/bland_altman_stats.csv`).

### 12.4 From many cells to an ECG

The ECG is a whole-heart signal, and its shape comes from differences between cells, not from any one cell. The QRS complex is the wave of depolarisation sweeping through the ventricles, carried by the His-Purkinje system. The T wave is the spread of repolarisation times across the wall and from apex to base. If every cell repolarised at the same moment the T wave would vanish. Mid-myocardial (M) cells have the longest action potentials, so the end of the T wave is set largely by the last tissue to recover.

The QT interval is therefore activation time plus the repolarisation time of the latest region. Two consequences matter here. First, electrical coupling between neighbouring cells blunts the differences. Real tissue shows about 30 to 40 ms of dispersion where uncoupled cells would show far more. In our uncoupled simulations the QT surrogate is the M-cell APD90 in every condition where that cell repolarises, with an implied dispersion of 101.9 ms. That is why the common "weighted average APD90" mapping has no mechanism behind it. Second, a dish monolayer has no activation sequence, no Purkinje system, no gradient across a wall and no autonomic input, so it has no analogue of the T wave at all.

## 13. The four corrections, and what each one does

### 13.1 Rate correction

Repolarisation takes longer when the heart beats slower, so a raw QT or FPD cannot be interpreted without the rate. The standard repair divides by a power of the cycle length: Ic = I / CL^α, with α = 1/3 for Fridericia and 1/2 for Bazett. We fit α on drug-free records only, baseline and placebo, because fitting on treated records would absorb part of the drug effect into the correction. The fitted exponents were 0.301 and 0.299, close to Fridericia's 1/3. That is the expected answer in healthy volunteers and a check on the code.

The dish needs this correction more than the clinic does. hiPSC-CMs beat spontaneously at roughly 29 to 37 beats per minute, and a drug that changes the beat rate changes the raw FPD for reasons that have nothing to do with repolarisation.

### 13.2 Hysteresis

QT does not follow a change in heart rate at once. It lags with a time constant of roughly one to two minutes, so the QT in any beat reflects a weighted average of the preceding minute or two of RR intervals and not that beat's own RR. The standard correction builds exactly that: an exponentially weighted average of preceding RR intervals, used in place of the instantaneous RR.

We could not do this, and the reason is built into the data. Each public record is a 10-second strip. That is 8 to 17 % of one adaptation time constant, so most of the RR history that determines the QT in the strip happened before the strip began. No weighting of 10 seconds of data can recover information that is not in it. The cost is measurable. Within a single window, QTcF computed on the slowest and the fastest beat differs by a median of 12.6 ms, and by more than 10 ms in 66.8 % of 238 records.

### 13.3 Autonomic tone

There is no standard "autonomic correction", and it helps to see why. Sympathetic and parasympathetic activity change both the heart rate and the repolarisation itself, and the second part is not removed by rate correction. Trials usually control it through design: supine rest before recording, matched times of day, and a placebo arm that is subtracted so that circadian and postural effects cancel.

The size of what is being controlled is not small. Ibutilide in healthy subjects with intact autonomic tone prolonged QTc by 43 ± 8 ms. After double autonomic blockade the same drug produced a further 68 ± 5 ms (Smith 2007). Tone attenuates drug-induced prolongation. A dish has no autonomic input, so it sits permanently in the state where a hERG blocker's effect is largest. Neither side can correct that away. It runs in the right direction to explain part of the dish's over-reading.

### 13.4 Irregular beats and the median beat

QT is only defined on a sinus beat that follows a reasonably regular preceding beat. An ectopic beat, or the compensatory pause after it, changes repolarisation for reasons unrelated to the drug. The usual handling is to exclude non-sinus beats and to require a stable preceding RR.

The median beat is the standard tool. Every beat in the window is aligned on a fiducial point, usually the R peak, and the median is taken across beats, sample by sample. It suppresses noise and is robust to a single odd beat, which is why it is preferred to the mean. One thing is often misstated. The median beat does not give the √N gain in signal-to-noise that averaging N beats would. That belongs to the mean. For Gaussian noise the median is about 1.25 times less efficient than the mean, and it is used for robustness and not for efficiency.

In the dish the analogous problem is harsher and cannot be averaged away. A well that develops an early afterdepolarisation has no measurable field potential duration at all, so the observation is lost and not merely noisy. It is lost exactly when the drug effect is largest. That is the informative missingness that rises from 0.9 % to 60.8 % across the concentration range.

## 14. Why one solution does not fit all three arms

The cascade is applied in the same form to all three arms, which is the point of the package. But the same form does not mean the same difficulty. Every stage means something different in each arm, and several problems have a fix in one arm that is meaningless or impossible in another.

| Stage | In vitro | In vivo | In silico |
|---|---|---|---|
| Rate | The preparation sets its own rate, 29 to 37 bpm, and drugs change it. The correction is large and unavoidable. | Rate is set by autonomic tone. The correction is modest at resting rates. | The pacing protocol imposes the rate. There is nothing to correct. |
| Baseline | The well before dosing. Available in principle, but not shipped in the public workbook. | The subject's own pre-dose record, within period. Available. | The drug-free run of the same model. Exact and noise-free. |
| Control | Vehicle wells, usually DMSO. | A time-matched placebo arm. | None needed, and none exists. |
| Hysteresis | Not applicable: the preparation is at steady state at its own rate. | Applicable and important, but unidentifiable in 10-second records. | Not applicable: pre-paced to steady state. |
| Autonomic tone | Absent entirely. It can be acknowledged but not corrected. | Present, and controlled by design. | Absent by construction. This is a modelling choice, not a limitation. |
| Missing data | Informative and severe: arrested wells, 0.9 % to 60.8 % with concentration. | Rare and close to random. | Informative but deterministic: from about 3.8× the mid-myocardial cell stops repolarising. |
| Exposure | Nominal concentration, reduced by an unmeasured amount through protein binding, plastic and membrane. | Measured free plasma concentration. | Exactly the concentration prescribed, with no uncertainty. |
| Dominant variance | Between site × platform × cell line. | Between subjects. | Parameter uncertainty, and only for hERG. |
| Endpoint | FPD, or APD90 in the optical arm, pooled in one column. | QT, J-Tpeak, Tpeak-Tend. | APD90 and qNet. |

Four consequences follow, and each one blocks a repair that works elsewhere.

Averaging fixes noise in the clinic and loses data in the dish. A noisy clinical beat can be averaged away with a median beat. An arrhythmic well produces no number at all, so there is nothing to average, and the observation is destroyed where the drug effect is largest.

Placebo subtraction removes common drift in the clinic and has no referent in a simulation. The model has no circadian rhythm, no meal and no change of posture. Applying the correction there would be meaningless and not merely conservative.

Hysteresis is a problem for the clinic alone. The dish and the model are at steady state, and an effect-compartment model would improve neither.

Exposure uncertainty runs in opposite directions. In silico the concentration is known exactly and the response is uncertain. In the clinic the concentration is measured with some error but the response is observed. In the dish the concentration is the least certain quantity in the whole chain, and the entire comparison rests on it. That is why the shared axis is the weakest link, and why it gets its own paragraph in section 7.

In practice only the concentration axis can be shared. Everything else has to be handled arm by arm and then compared. A single harmonisation procedure applied to all three would be wrong even if it were convenient.

## 15. Do the three systems agree?

Asked in that bare form, the question is not well posed. Asked precisely, it has a standard answer, and the precise version is the one the field uses.

Start with why the bare question fails. The three arms measure different quantities, so agreement in milliseconds is not defined. Agreement also means nothing without an exposure attached. In an FDA database of 150 drug applications, the hERG, action potential duration and in vivo QTc assays were specific but not sensitive (Park 2018). Sensitivity was low (0 to 0.33) at 1 to 10 times the free clinical exposure. It rose to 0.67 to 0.80 at higher multiples for hERG and in vivo QTc, while the action potential assay stayed at or below 0.53. Specificity was high throughout, reaching 0.89 to 1.0 at multiples up to 30. The same assay can therefore look concordant or discordant depending only on where it is evaluated. (The action potential assay in that study is the conventional non-clinical one, not the hiPSC-CM assay.) And the reference standard is itself a surrogate: a thorough-QT outcome measured in healthy volunteers, none of whom had an arrhythmia.

What the field measures instead is diagnostic accuracy against a named clinical endpoint: sensitivity, specificity, predictive values and likelihood ratios. Gintant (2022) reviews these measures, including the gain from combining two independent assays. The ICH E14/S7B Questions and Answers build on the same idea with a "double negative", meaning a negative hERG assay and a negative in vivo QTc study, which can reduce the number of thorough QT studies. Strauss and colleagues (2021) reviewed an industry proposal (Vargas et al.) that double-negative drugs might not need detailed clinical QTc evaluation. They did not accept its exact probability figures, because of limits in the underlying datasets, but showed that even modest assay performance can still leave a low probability of torsades after a double negative. They also note that the 20 % pre-test probability of QTc prolongation assumed in that proposal is consistent with FDA data in which 19 % of drugs prolonged QTc. The rule depends on two negative results agreeing. It says nothing about the systems producing the same numbers.

### 15.1 Would a relative scale fix it?

Partly, and the part that fails is informative. We tested it directly (`results/relative_scale_test.csv`). On a percent-of-baseline scale the dish-to-clinic gap narrows from 1.81× to 1.48× but does not close. Normalised instead to a reference drug measured in the same system, the two pure hERG blockers agree almost exactly. Dofetilide-normalised quinidine is 1.023 in the dish against 1.025 in the clinic, a ratio of 0.999. Ranolazine is still 1.92× off and verapamil inverts. A relative scale removes the units problem and the common multiplicative factor. It leaves the drug-specific disagreement alone, and no normalisation will remove that residue.

Two cautions apply. The dish baseline that a percent scale needs is not in the public workbook, so we imported one from an independent study (the median control FPDc of the Millard 2018 studies, 482 ms). And verapamil's clinical interval spans zero, so any ratio with it in the denominator has unbounded variance.

### 15.2 Would a Bland-Altman plot settle it?

Only where the two methods measure the same quantity. For dish against clinic it is not appropriate. The average of an FPD and a QT estimates nothing, and four drugs cannot support limits of agreement. Where it is appropriate, it is informative (Figure 7). Optical APD90 against electrode FPD, over 97 paired drug and concentration cells, gives a bias of +34.7 ms and limits of agreement from −227 to +297 ms. There is strong proportional bias (slope 0.800, p < 0.001; 0.502 with the single extreme cell removed). A band 524 ms wide, against drug effects of 10 to 130 ms, is the sharpest test we have of the source study's description of the two as equivalent measures. One caution: the two dish methods were never applied to the same wells, so this measures the total disagreement between the two arms and not the endpoint effect alone. Our tangent T-end against the study's own, on 220 records in 22 subjects, gives a bias of −9.1 ms with limits of −64 to +46 ms and no proportional bias (p = 0.49).

![Figure 7. Bland-Altman agreement. (a) Optical APD90 against electrode FPD over 97 paired cells, with a wide band and strong proportional bias. (b) Our tangent T-end against the study's own on 22 subjects, with a narrow band and no proportional bias. Limits of agreement use between- and within-cluster variance components, because several concentrations share a drug and several records share a subject.](../results/bland_altman.png)

Limits of agreement use the closed-form unbalanced one-way variance components, not a naive SD, because treating clustered observations as independent would understate the spread. A random-intercept model was tried first and did not converge on the clustered optical comparison.

### 15.3 The useful form of the question

The productive question is what each system adds that the others do not. Three systems that always agreed would make two of them redundant. The regulatory framework values the combination because the assays fail in different directions.

## 16. What the data cannot support, and what we would want

### 16.1 Random effects that cannot be identified

Below the dataset level in vitro, the workbook gives one ΔΔFPDc per well and no plate identifier, well position, electrode count, timestamp, cell lot or passage number. Plate-to-plate variation and edge effects, both well documented in plate assays, have nowhere to go in the model. The finest grouping the data supports is dataset = site × cell line × platform, which is what we fit, and everything below it falls into the residual. The dish intervals are 25 to 106 ms wide and we attribute that to variation between datasets, but we cannot exclude plate structure as part of the cause.

Across channels in silico, `data/hill28.csv` carries IC50 and Hill coefficient as single point estimates for all seven currents. Only hERG has a distribution (`data/herg_boot28.csv`, 20 bootstrap sets per drug drawn from the 2,000 published). Every qNet interval in this work therefore carries hERG binding uncertainty alone and treats the other six channels as exactly known. The true uncertainty on qNet is larger than we report, by an amount the public data cannot show.

### 16.2 The data we would ask for, in order of value

1. Raw, uncorrected well-level FPD and beat period, with plate and well coordinates and read timestamps. This one item would allow cascade validation on real wells, the plate random effect, percent-change normalisation and the equilibration-drift question.
2. Measured free concentration per well, per drug and per site.
3. A paired endpoint set: the same cells read by both platforms in the same medium.
4. A reason code for every missing well, separating arrhythmia from signal loss, from never-beating, from quality exclusion.
5. IC50 and Hill uncertainty for all seven channels, not hERG alone.
6. Clinical recordings longer than 10 seconds that carry the preceding RR history, so hysteresis becomes identifiable.
7. Beat-level clinical annotations instead of one median beat per record.
8. More open single-agent concentration-QTc datasets. Four drugs is the whole bridge.
9. A subject-matched dish-and-clinic study larger than 16 people.
10. Any torsades outcomes at all. The corpus has none.

## 17. Literature summary

This section says what each key paper contributes, what we used it for, and the one caution that matters when leaning on it. We read most of these papers in full. Where we relied only on an abstract or an excerpt, the entry says so.

### 17.1 The in-vitro evidence

Blinova and colleagues (2018, Cell Reports) report the CiPA myocyte study: 28 drugs at ten sites, two commercial hiPSC-CM lines and six platforms. It is the largest public dataset of its kind and the whole dish side of this work. It ships ΔΔFPDc only, already corrected for baseline, vehicle and rate, and it pools electrode FPD with optical APD90 in one column while calling them equivalent measures.

Blinova and colleagues (2016, Toxicological Sciences) is the earlier pilot, which showed that the assay separates torsadogenic from non-torsadogenic compounds. We use it only as precedent for drug-level use.

Blinova and colleagues (2019, Clinical and Translational Science) is the "clinical trial in a dish". Cells were reprogrammed from 16 of the volunteers who took the drugs, so the dish and the clinic refer to the same people. It found no within-person correlation (ρ = −0.11, p = 0.70 for dofetilide). Those subjects' clinical responses to two drugs did not correlate with each other either (ρ = 0.10), so part of the null is on the clinical side. This is why the package stays at drug level.

Millard and colleagues (2018, Toxicological Sciences) ran a blinded cross-site MEA pilot. Its supplementary table gives baseline beat period, FPD and corrected FPD for 18 studies, and it is the only public dataset on which our rate correction could be checked against an independent implementation (17 of 18 within 5 ms, median 0.65 ms). Study 9 is off by 22.0 ms and we cannot explain it.

Schocken and colleagues (2018) studied media effects. Serum-containing against serum-free medium changes baseline FPDc by 104 ms (575 ± 59 against 471 ± 44) and slows beating by a third. This is why the confound at the optical site cannot even be given a sign. It also overturned one of our claims: of 18 drugs that induced arrhythmias, 8 did so at lower concentrations with serum present.

Guerrelli and colleagues (2024) measured temporal drift. Beat rate rose 22.6 % and FPD fell 7.7 % over a 20-minute equilibration, which is the same order as the drug effects we are measuring. Ernault and colleagues (2024) explain what a field potential is and why the criteria for reading it differ from those for an electrogram. We rely on that paper for section 12.3.

### 17.2 The clinical evidence

Johannesen and colleagues (2014) report the ECGRDVQ trial: dofetilide, quinidine, ranolazine, verapamil and placebo in 22 subjects, crossover, with concentrations. Three of our four bridge drugs come from it. Johannesen and colleagues (2015) report ECGDMMLD, on late-sodium block. We reproduce its dofetilide-alone result to within 0.2 ms at two time points. That resolved what had looked like our own 6 to 7 ms error, which turned out to be a disagreement between two published analyses of the same trial.

Vicente and colleagues (2015) analysed T-wave morphology in ECGRDVQ, and its Table 2 is our main validation target (8 of 8). Vicente and colleagues (2016) and (2019) cover the late-sodium biomarker work and the CiPA phase I validation study. The 2019 paper is the source of the ECGCIPA data and of the prespecified J-Tpeak power-law exponent of 0.58, which we confirmed on 5,740 records. Its Table 3 and the primary paper disagree on dofetilide, and we match the primary paper. Johannesen and colleagues (2016) describe the automated J-Tpeak and Tpeak-Tend algorithm behind the interval tables we use.

Hnatkova and colleagues (2019) give the linear J-Tpeak rate correction, the alternative to the power law. Per record the two differ by a median of 0.44 ms, and at drug level by 1.2 ms or less. Malik and colleagues (2019) explain why population-wide rate corrections are inadequate and individual QT/RR profiles matter, which is the basis for our per-subject exponent option. Garnett and colleagues (2016) describe the operating characteristics of linear concentration-QT models, the design our clinical fits follow.

Strodthoff and colleagues (2023) give PTB-XL+, where three algorithms measured the same 21,795 ECGs. It is our outside evidence that interval extraction carries a rate-dependent bias of about 6 ms across a 300 ms change in RR. Smith and colleagues (2007) showed that autonomic tone attenuates drug-induced QT prolongation: ibutilide prolonged QTc by 43 ± 8 ms with tone intact and by a further 68 ± 5 ms under double blockade. That is the quantitative basis for expecting a denervated preparation to over-read.

### 17.3 The in-silico and ion-channel evidence

O'Hara and colleagues (2011) is the human ventricular model behind our 408 simulations. They report an endocardial APD of about 275 ms in human preparations at a 1,000 ms cycle length, matched by the model at 273 ms. Our control endocardial APD90 is 265.5 ms, which is 2.7 % below. The text does not state the epicardial and M-cell control values numerically, so those two remain unchecked.

Dutta and colleagues (2017) reparameterised the model for CiPA and defined qNet. Li and colleagues (2017) combined hERG binding kinetics with multichannel pharmacology, which is the basis of the trapping formulation, and Li and colleagues (2019) assessed the in-silico metric formally under CiPA. Crumb and colleagues (2016) tested 30 drugs against the CiPA channel panel, and its Table 1 is our source of free Cmax. It does not tabulate protein binding.

Alvarez Baron and colleagues (2025) describe patch-clamp practice, including concentration verification in the recording rig. Measured losses ranged from 6 % (ranolazine) to 58 % (vanoxerine), with verapamil at 41 % and chloroquine at 19 %. This is the quantitative anchor for section 7, with the caution that it comes from perfusion rigs and not from static MEA wells. Pritchard and colleagues (1985) measured bepridil binding to plasma protein by equilibrium dialysis and showed that it saturates.

### 17.4 Framework and translation

Sager and colleagues (2014), Colatsky and colleagues (2016) and Fermini and colleagues (2016) describe the CiPA paradigm, its rationale and its progress. They give context and no data. Park and colleagues (2018, British Journal of Pharmacology) analysed 150 drug applications in an FDA database. The paper gives the sensitivity and specificity figures quoted in section 15 and found that the in vivo QTc assay gave the most robust predictive values, while the conventional action potential assay added little. We read this one in full.

Gintant (2022) is a single-author review of how to judge whether a non-clinical assay translates. We read only the abstract. Strauss and colleagues (2021) review the ICH E14/S7B updates and the double-negative proposal. We read it in full.

Stillitano and colleagues (2017) and Shinozawa and colleagues (2017) report positive subject-level correlations between hiPSC-CM responses and clinical ones, and Blinova 2019's null has to be read against them. We know them only through Blinova 2019. The literature is unsettled and this package does not try to settle it. Redfern and colleagues (2003) proposed the provisional safety-margin framework. The paper is closed access, and one excerpt read through a search result is the basis for our terfenadine exposure note.

## 18. Were the processing choices right?

Every number here depends on four choices that could have gone another way: the rate-correction exponent, the shape of the concentration predictor, the random-effect structure, and whether the placebo adjustment is formed in advance or estimated inside the model. This section writes each model out, refits the headline quantity under the alternative, and reports what moves. The code is `scripts/methodology_sensitivity.py` and the output is `results/methodology_sensitivity.csv`.

### 18.1 The models, written out

In vivo, as implemented here. For subject i and post-dose record j, with concentration C centred at the geometric-mean observed Cmax, C*:

    ΔΔQTcF_ij = β0 + β1 (C_ij − C*) + b_i + ε_ij,   b_i ~ N(0, σ²_b),  ε_ij ~ N(0, σ²)

Here β0 is the effect at C*. It is reported with a two-sided 90 % interval that uses a t quantile on (number of subjects − 1). The double delta is formed before the model: baseline subtracted within subject and period, then the time-matched placebo arm subtracted.

In vivo, as the field prescribes. The prespecified linear mixed-effects model in the concentration-QTc white paper (Garnett 2018) uses the single delta as the response and estimates the placebo adjustment inside the model:

    ΔQTcF_ij = θ0 + θ1 C_ij + θ2 (QTc0_i − mean QTc0) + θ3 TRT_ij
               + Σ_k θ4k 1{time = k} + b_i + b'_i C_ij + ε_ij

Placebo records enter at C = 0 and TRT = 0, and the reported quantity is ΔΔQTc(C) = θ1 C + θ3. Treatment and post-baseline time point are fixed effects, the subject has a random intercept and a random slope, and the centred baseline QTc enters as a covariate.

In vitro. For dataset d (site × cell line × platform, 15 of them) and well w, with m the multiple of free Cmax:

    ΔΔFPDc_dw = γ0 + γ1 (log10 m_dw − log10 m*) + u_d + ε_dw,   u_d ~ N(0, τ²)

The random effect is the dataset and not the well, because nothing below the dataset level can be identified (section 16.1).

In silico. There is no mixed model here, and that is not an omission. A simulation has no sampling error: a given parameter vector returns one answer. The uncertainty that exists is parametric, y_b = g(θ_b) over bootstrap draws θ_b of the hERG binding parameters, and it covers one channel of seven. The band in Figure 8c is the range across the three cell types. It is a spread and not a confidence interval, and the figure says so.

![Figure 8. Dofetilide concentration-response in each arm with 90 % bands, read at the geometric-mean observed free Cmax and annotated with the interval there. (a) In vitro: +132 ms [93, 171]. (b) In vivo: +73 ms [69, 78] against the ICH E14 10 ms line. (c) In silico: +50 ms [41, 50], where the band is the range over endocardial, epicardial and mid-myocardial cells and not a confidence interval.](../results/confidence_bands.png)

### 18.2 Does one cascade make sense?

In form, yes. In parameters, no. The same ordered sequence of rate, baseline, control, model and shared scale is what makes the comparison auditable, and enforcing it in code stops an endpoint being used before it has been corrected. But every parameter inside the sequence has to be specific to the arm, and two of the five steps are empty in at least one arm (section 14). A single cascade with a single set of parameters would be wrong.

### 18.3 What if the power law differs by arm?

It does, and this is where the largest mistake could have been made. Clinically, Fridericia is defensible on our own data, since the exponents fitted on drug-free records are 0.301 and 0.299. The choice still matters:

| Drug | Fridericia (used) | Bazett | Individual exponent |
|---|---|---|---|
| Dofetilide | 73.1 [68.5, 77.7] | 72.4 [67.3, 77.5] | 72.3 [67.4, 77.3] |
| Quinidine | 74.9 [69.3, 80.5] | 85.4 [78.4, 92.4] | 78.6 [72.1, 85.1] |
| Ranolazine | 11.2 [7.6, 14.8] | 14.4 [9.6, 19.2] | 12.4 [7.8, 17.0] |
| Verapamil | 3.5 [−0.1, 7.0] | 10.5 [6.0, 15.0] | 5.8 [1.7, 9.8] |

Bazett flips the ICH E14 decision for verapamil. Under Fridericia the upper bound is 7.0 ms and the drug is negative. Under Bazett it is 15.0 ms and the drug is positive. This is the textbook failure of Bazett, which over-corrects when heart rate moves, and verapamil slows the heart. It is still instructive to see it change a regulatory call on real data. The result is asserted in the test suite.

In vitro we could not make this choice, because the workbook arrives already corrected. Fridericia is the usual standard for field potentials. The Axion induced-model standard defines FPDc as FPD divided by beat period to the power 0.33, a convention adopted to reduce the influence of a wide range of interspike intervals. But the standard is not unchallenged. A PLOS ONE characterisation study of hiPSC-CM electrophysiology (doi:10.1371/journal.pone.0167348), which we read only as a search excerpt, notes that there is no evidence that clinical correction formulae apply to FPD from hiPSC-CMs of varying beat rate. It finds Fridericia and Bazett both fairly effective, with Fridericia likely under-correcting and Bazett over-correcting at slow interspike intervals. A separate optimisation study, also seen only as an abstract, developed and qualified correction algorithms by moving the beat rate on purpose with the If inhibitor ZD7288 and isoproterenol.

So the best method in vitro is probably not to assume an exponent but to qualify one for each preparation with a rate-changing tool. We could not do that. The only check available was indirect. The exponent implied by the Millard 2018 baseline data is 0.324 against Fridericia's 0.333, which is reassuring but comes from a different study's cells.

### 18.4 What if the predictor shape were the same in both arms?

This choice has the largest consequence, and it cannot be shared.

| Drug | log10 multiples (used) | linear in multiples |
|---|---|---|
| Dofetilide | +132.0 [93.1, 170.8] | +96.2 [60.9, 131.6] |
| Quinidine | +121.7 [76.6, 166.7] | +126.8 [85.5, 168.1] |
| Ranolazine | +39.5 [26.4, 52.6] | +33.8 [21.2, 46.4] |
| Verapamil | −87.2 [−99.3, −75.1] | −43.4 [−53.8, −33.0] |

Forcing a linear predictor halves verapamil's effect and cuts dofetilide's by about a quarter. Tested in-vitro concentrations span up to 3,500-fold, and a linear fit is dominated by its top point. Forcing a log predictor on the clinical side would be wrong in the other direction. Clinical concentrations span less than two decades, and the ICH E14 decision is defined on a linear concentration scale. Each arm needs a different predictor, for reasons internal to that arm.

### 18.5 What if we pooled differently?

Across datasets in vitro, treating the dataset as a fixed effect and not a random one moves dofetilide from +132.0 to +118.3 ms, quinidine from +121.7 to +123.6, ranolazine from +39.5 to +39.3 and verapamil from −87.2 to −86.2. The random-effect version shrinks datasets toward the common mean and supports a statement about sites in general. The fixed-effect version describes these 15 datasets and nothing else. We want to generalise, so the random effect is the right default, but the 13.7 ms shift for dofetilide shows the choice has a cost.

Across studies, we deliberately do not pool ECGRDVQ with ECGCIPA. At matched concentration they differ by 6 to 7 ms for the same drugs, and we cannot explain it (section 4.7). Pooling would average that discrepancy into a number that looks more precise than the evidence allows.

### 18.6 Where our clinical model departs from the standard

There are four departures.

1. We form the double delta in advance instead of estimating the placebo adjustment inside the model through treatment and time effects. We tested this directly. The white-paper-style model gives 69.8 against our 73.1 (dofetilide), 71.7 against 74.9 (quinidine), 11.7 against 11.2 (ranolazine) and 3.5 against 3.5 (verapamil). The largest gap is 3.4 ms and no decision changes.
2. We leave out the centred baseline QTc covariate that the white paper requires. This is a real gap and not a tested and rejected choice. It is the departure we would fix first.
3. We use a random intercept only. Adding a random slope moves the effects by 0.0 to 3.0 ms and changes no decision. An earlier test over nine published comparisons moved the mean absolute gap from 3.16 to 3.04 ms.
4. We use no categorical time effect. That is partly defensible, since time-matched baseline adjustment is argued to reduce diurnal variation enough that categorical time effects are not needed. We did not test it.

If starting again, the best method is the prespecified model as written: ΔQTcF as the response, with concentration, centred baseline, treatment and time as fixed effects, a random intercept and slope for each subject, and the double delta taken from model contrasts and not formed in advance. We now know what the shortcut cost on these four drugs, at most 3.4 ms and no change to any conclusion. That is the honest way to report a simplification: by measuring how much it matters and not by asserting that it does not.

### 18.7 Is the concentration-response linear at all?

The in-silico arm fits nothing. This is worth saying because Figure 2 looks like a steeply curving dose-response and invites the question. The simulation returns one action potential per parameter vector. There is no sampling, no residual and no regression. The in-silico effect we quote at 1.07× free Cmax was simulated directly, so not even interpolation is involved.

The curve in that figure is also partly an effect of the log axis. Plotted against concentration and not log concentration, the simulated response is close to a straight line across the exposures that matter. A linear fit of dAPD90 on multiples of free Cmax gives R² = 0.9993 for dofetilide (largest residual 2.2 ms), 0.9894 for quinidine, 0.9820 for ranolazine and 0.9089 for verapamil. The convexity is real but supratherapeutic. Extending dofetilide to 30× in the endocardial cell, which does keep repolarising, drops the linear fit to R² = 0.9387 with a largest residual of 70 ms. So the model is nearly linear through the therapeutic range and accelerates well above it, which is also where it eventually fails to repolarise.

Where linearity is assumed, in the clinical model, we tested it. Linear is the field default, but the same guidance requires the assumption to be checked and an Emax model fitted when it fails. Study protocols routinely prespecify that fallback. We fitted linear, log-linear and Emax forms with the same random intercept, by maximum likelihood so the fixed structures can be compared by BIC (`results/model_shape.csv`).

| Drug | Best by BIC | Linear ΔBIC | Effect at Cmax, linear | Effect at Cmax, best | E14 call changes? |
|---|---|---|---|---|---|
| Dofetilide | linear | 0.0 | +73.1 [68.6, 77.6] | +73.1 | no |
| Quinidine | log-linear | 19.3 | +74.9 [69.4, 80.4] | +73.2 [68.0, 78.4] | no |
| Ranolazine | log-linear (tie) | 0.4 | +11.2 [7.7, 14.8] | +11.6 [8.1, 15.2] | no |
| Verapamil | Emax | 7.1 | +3.5 [0.0, 6.9] | +4.0 [0.6, 7.4] | no |

Linear is the best form for only one of the four drugs. Quinidine clearly prefers a curved form, and a recent paper on concentration-QT analysis beyond the linear model (Cellière and colleagues, PMC12065759, which we read only as a search excerpt) reports the same for quinidine, with Emax and log-linear fitting better than linear by BIC. Verapamil prefers Emax. It is also the drug whose two misses against the published analysis we had already put down to fitting a linear model to a curved response.

The consequence is small. The effect at the geometric-mean Cmax moves by at most 6 ms across the three shapes, and no ICH E14 decision changes for any of the four drugs. Compare that with the rate-correction choice in section 18.3, where swapping Fridericia for Bazett flips verapamil's decision outright. Curvature is the more interesting modelling question and the less consequential one.

If curvature mattered, we would try Emax first, f(C) = Emax·C/(C + EC50), which is what protocols prespecify. Sigmoid Emax adds a Hill coefficient and is what the receptor pharmacology implies. Log-linear, f(C) = p1·log(1 + C/p2), sits in between. In vitro we already use a nonlinear transform, since the predictor is log10 of the multiple, which is a local approximation to the middle of a sigmoid over a 3,500-fold span.

The same paper reports that adding hysteresis to the linear model lowered BIC by 105.8 for ΔΔQTc and by 109.6 for ΔQTc. We saw only a search excerpt of that table. We cannot say which drugs it covers, and the excerpt does not show the Emax-versus-linear gap, so we do not compare the two sizes. In our own fits the shape of the curve moves the effect at Cmax by at most 6 ms. Hysteresis is exactly what the 10-second public records cannot support (section 13.2). Our judgement, and not the paper's finding, is that if one modelling improvement were available to this project, it would be longer recordings and not Emax.

## 19. Availability

The code is the `dish2ecg` package (142 tests, version 0.2.2) under the MIT licence, copyright Rejaul Karim. The repository is at https://github.com/rejaul/dish2ecg. `results/README.md` maps every result file to the script that produced it and gives the order for regenerating them. The well-level dish table is not redistributed. `scripts/build_wells.py` rebuilds it from the original workbook.

## Glossary

- QT interval: the time on an ECG from the start of ventricular depolarisation to the end of the T wave. A long QT predisposes to torsades de pointes.
- QTc, QTcF: QT corrected for heart rate. The F stands for Fridericia's correction, which divides by RR to the power 1/3.
- FPD: field potential duration, the dish equivalent, measured with electrodes under the cells.
- APD90: action potential duration to 90 % repolarisation, measured in single cells or with optical probes.
- Δ and ΔΔ: the change from baseline, and the change from baseline minus the placebo or vehicle change.
- EAD: early afterdepolarisation, an abnormal extra depolarisation during repolarisation. In the dish it marks arrhythmic behaviour.
- Free Cmax: the peak concentration of unbound drug in plasma at a clinical dose.
- Nominal concentration: in a dish experiment, the concentration worked out from the amount of drug added to the well. The concentration the cells actually see is lower, because drug binds medium protein, sticks to plastic and tubing, and partitions into the cells.
- Fraction unbound (fu): the share of drug in plasma that is not bound to protein.
- CiPA: the Comprehensive in vitro Proarrhythmia Assay, a regulatory initiative that proposed combining ion-channel data, a computer model and hiPSC-CM data.
- qNet: the CiPA in-silico score, the net charge carried by the main ion currents during a beat, in µC/µF. Lower values mean higher torsades risk.
- Hysteresis: the delay with which QT follows a change in heart rate.
- AUC: area under the ROC curve. 0.5 is chance and 1.0 is perfect separation.

## References

- Alvarez Baron C, et al. (2025) Reassessing the conduct of patch clamp cardiac ion channel assays to improve nonclinical data translation to clinical ECG changes and proarrythmia risk. *Clin Transl Sci* 18:e70451. doi:10.1111/cts.70451
- Axion BioSystems. AIMS #CM01, iPSC model standard. https://www.axionbiosystems.com/axion-ipsc-model-standards-aims/aims-cm01
- Blinova K, et al. (2016) Comprehensive translational assessment of human-induced pluripotent stem cell derived cardiomyocytes for evaluating drug-induced arrhythmias. *Toxicol Sci* 155:234–247. doi:10.1093/toxsci/kfw200
- Blinova K, et al. (2018) International multisite study of human-induced pluripotent stem cell-derived cardiomyocytes for drug proarrhythmic potential assessment. *Cell Reports* 24:3582–3592. doi:10.1016/j.celrep.2018.08.079
- Blinova K, et al. (2019) Clinical trial in a dish: personalized stem cell-derived cardiomyocyte assay compared with clinical trial results for two QT-prolonging drugs. *Clin Transl Sci* 12:687–697. doi:10.1111/cts.12674
- Cellière G, Krause A, Bonnefois G, Chauvin J. Beyond the linear model in concentration-QT analysis. PMID 40347307; PMC12065759 (journal and year not confirmed; read as a search excerpt only)
- Colatsky T, et al. (2016) The Comprehensive in Vitro Proarrhythmia Assay (CiPA) initiative: update on progress. *J Pharmacol Toxicol Methods* 81:15–20. doi:10.1016/j.vascn.2016.06.002
- Crumb WJ, et al. (2016) An evaluation of 30 clinical drugs against the comprehensive in vitro proarrhythmia assay (CiPA) proposed ion channel panel. *J Pharmacol Toxicol Methods* 81:251–262. doi:10.1016/j.vascn.2016.03.009
- Dutta S, et al. (2017) Optimization of an in silico cardiac cell model for proarrhythmia risk assessment. *Front Physiol* 8:616. doi:10.3389/fphys.2017.00616
- Electrophysiological characteristics of human iPSC-derived cardiomyocytes for the assessment of drug-induced proarrhythmic potential. *PLoS ONE*, e0167348. doi:10.1371/journal.pone.0167348 (first author not verified; read as a search excerpt only)
- Ernault AC, et al. (2024) Interpretation of field and LEAP potentials recorded from cardiomyocyte monolayers. *Am J Physiol Heart Circ Physiol* 326:H800–H811. doi:10.1152/ajpheart.00463.2023
- Fermini B, et al. (2016) A new perspective in the field of cardiac safety testing through the comprehensive in vitro proarrhythmia assay paradigm. *SLAS Discov* 21:1–11. doi:10.1177/1087057115594589
- Garnett C, et al. (2016) Operational characteristics of linear concentration-QT models. *Clin Pharmacol Ther* 100:170–178. doi:10.1002/cpt.361
- Garnett C, et al. (2018) Scientific white paper on concentration-QTc modeling. *J Pharmacokinet Pharmacodyn* 45:383–397. doi:10.1007/s10928-017-9558-5
- Gintant G (2022) Assessing the fidelity of translation of non-clinical assays: a Pharma perspective. *Br J Pharmacol*. doi:10.1111/bph.15796 (abstract only)
- Goldberger AL, et al. (2000) PhysioBank, PhysioToolkit, and PhysioNet. *Circulation* 101:e215. doi:10.1161/01.cir.101.23.e215
- Guerrelli D, et al. (2024) hiPSC-CM electrophysiology: impact of temporal changes and study parameters on experimental reproducibility. *Am J Physiol Heart Circ Physiol* 327:H12–H27. doi:10.1152/ajpheart.00631.2023
- Hnatkova K, et al. (2019) Heart rate correction of the J-to-Tpeak interval. *Sci Rep* 9. doi:10.1038/s41598-019-51491-4
- Johannesen L, et al. (2016) Automated algorithm for J-Tpeak and Tpeak-Tend assessment of drug-induced proarrhythmia risk. *PLoS ONE* 11:e0166925. doi:10.1371/journal.pone.0166925
- Johannesen L, et al. (2014) Differentiating drug-induced multichannel block on the electrocardiogram: randomized study of dofetilide, quinidine, ranolazine, and verapamil. *Clin Pharmacol Ther* 96:549–558. doi:10.1038/clpt.2014.155
- Johannesen L, et al. (2015) Late sodium current block for drug-induced long QT syndrome: results from a prospective clinical trial. *Clin Pharmacol Ther* 99:214–223. doi:10.1002/cpt.205
- Li Z, et al. (2017) Improving the in silico assessment of proarrhythmia risk by combining hERG channel-drug binding kinetics and multichannel pharmacology. *Circ Arrhythm Electrophysiol* 10. doi:10.1161/CIRCEP.116.004628
- Li Z, et al. (2019) Assessment of an in silico mechanistic model for proarrhythmia risk prediction under the CiPA initiative. *Clin Pharmacol Ther* 105:466–475. doi:10.1002/cpt.1184
- Malik M, et al. (2019) Implications of individual QT/RR profiles, part 1. *Drug Saf* 42:401–414. doi:10.1007/s40264-018-0736-1
- Millard D, et al. (2018) Cross-site reliability of human induced pluripotent stem cell-derived cardiomyocyte based safety assays using microelectrode arrays: results from a blinded CiPA pilot study. *Toxicol Sci* 164:550–562. doi:10.1093/toxsci/kfy110
- O'Hara T, et al. (2011) Simulation of the undiseased human cardiac ventricular action potential: model formulation and experimental validation. *PLoS Comput Biol* 7:e1002061. doi:10.1371/journal.pcbi.1002061
- Optimizing rate correction of field potential duration, a biomarker for QT risk assessment, in human iPSC-cardiomyocytes. Abstract only (ResearchGate record 276097641); authors, venue and year not verified
- Park E, Gintant GA, Bi D, et al. (2018) Can non-clinical repolarization assays predict the results of clinical thorough QT studies? Results from a research consortium. *Br J Pharmacol*. doi:10.1111/bph.14101
- Pritchard JF, et al. (1985) Plasma protein binding of bepridil. *J Clin Pharmacol* 25:347–353. doi:10.1002/j.1552-4604.1985.tb02853.x
- Redfern WS, et al. (2003) Relationships between preclinical cardiac electrophysiology, clinical QT interval prolongation and torsade de pointes for a broad range of drugs: evidence for a provisional safety margin in drug development. *Cardiovasc Res* 58:32–45. doi:10.1016/S0008-6363(02)00846-5 (closed access; one excerpt read via a search result only)
- Sager PT, et al. (2014) Rechanneling the cardiac proarrhythmia safety paradigm: a meeting report from the Cardiac Safety Research Consortium. *Am Heart J* 167:292–300. doi:10.1016/j.ahj.2013.11.004
- Schocken D, et al. (2018) Comparative analysis of media effects on human induced pluripotent stem cell-derived cardiomyocytes in proarrhythmia risk assessment. *J Pharmacol Toxicol Methods* 90:39–47. doi:10.1016/j.vascn.2017.11.002
- Shinozawa T, et al. (2017) Recapitulation of clinical individual susceptibility to drug-induced QT prolongation in healthy subjects using iPSC-derived cardiomyocytes. *Stem Cell Reports* 8:226–234. doi:10.1016/j.stemcr.2016.12.014
- Smith AH, et al. (2007) Autonomic tone attenuates drug-induced QT prolongation. *J Cardiovasc Electrophysiol* 18:960–964. doi:10.1111/j.1540-8167.2007.00901.x
- Stillitano F, et al. (2017) Modeling susceptibility to drug-induced long QT with a panel of subject-specific induced pluripotent stem cells. *eLife* 6. doi:10.7554/elife.19406
- Strauss DG, Wu WW, Li Z, Koerner J, Garnett C (2021) Translational models and tools to reduce clinical trials and improve regulatory decision making for QTc and proarrhythmia risk (ICH E14/S7B updates). *Clin Pharmacol Ther*. doi:10.1002/cpt.2137
- Strodthoff N, et al. (2023) PTB-XL+, a comprehensive electrocardiographic feature dataset. *Sci Data* 10. doi:10.1038/s41597-023-02153-8
- Vicente J, et al. (2015) Comprehensive T wave morphology assessment in a randomized clinical study of dofetilide, quinidine, ranolazine, and verapamil. *J Am Heart Assoc* 4. doi:10.1161/JAHA.114.001615
- Vicente J, et al. (2016) Electrocardiographic biomarkers for detection of drug-induced late sodium current block. *PLoS ONE* 11:e0163619. doi:10.1371/journal.pone.0163619
- Vicente J, et al. (2019) Assessment of multi-ion channel block in a phase I randomized study design: results of the CiPA phase I ECG biomarker validation study. *Clin Pharmacol Ther* 105:943–953. doi:10.1002/cpt.1303
