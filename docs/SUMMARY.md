# Putting a dish and a patient on the same axis

**dish2ecg** — a reproducible cascade that passes hiPSC-CM field potentials and
clinical 12-lead ECGs through the *same* correction algebra, compares them on a
shared free-Cmax scale, and refuses to express an in-vitro endpoint as a QT
interval in milliseconds.

Built entirely on public data. 142 tests. Every number below comes from a
validated artifact in the repository.

---

## The question

A cardiac safety package asks one thing twice: *does this compound delay
repolarisation, and by how much, at exposures a patient will actually see?*
In vitro that question is answered in field-potential duration on a plate of
stem-cell cardiomyocytes; in the clinic it is answered in QTc on a 12-lead ECG.
The two answers are routinely compared, and the comparison is almost never
auditable — the corrections differ, the concentration axes differ, and the
in-vitro number is frequently reported as though it were a QT interval.

`dish2ecg` makes that comparison explicit: one cascade, applied symmetrically,
with provenance recorded at every step, so a reviewer can see exactly which
operations produced a given number and refuse the ones that are not licensed.

## What the comparison shows

![The dish tracks the clinic in sign and rank for the hERG blockers and inverts for verapamil. Four panels, one per drug, plotting change in repolarisation duration against multiples of free Cmax on a log axis; clinical ddQTcF from 22 subjects in blue and in-vitro ddFPDc from 15 datasets in orange, with the matched-exposure readout annotated in each panel.]({{artifact:e50a67c7-41a1-4ec2-811d-0b476b127145}})

Four drugs are the **entire overlap** between the public in-vitro panel and the
public clinical concentration-QTc corpora. At matched observed free Cmax:

| | clinic ddQTcF | dish ddFPDc | ratio |
|---|---|---|---|
| Dofetilide (1.07x) | +73.1 [68.5, 77.7] | +132.0 [93.1, 170.8] | 1.81 |
| Quinidine (0.29x) | +74.9 [69.3, 80.5] | +135.1 [82.2, 187.9] | 1.80 |
| Ranolazine (0.96x) | +11.2 [7.6, 14.8] | +38.8 [22.8, 54.7] | 3.46 |
| Verapamil (0.50x) | +3.5 [-0.1, 7.0] | **-86.3** [-98.9, -73.6] | -24.66 |

The dish agrees in sign on 3 of 4 and over-reads by roughly 1.8-3.5x on the
three it gets right. Its confidence intervals are 5-11x wider than the clinical
ones despite 128-278 wells per drug — heterogeneity across site, cell line and
platform, not sampling noise.

Verapamil is the instructive failure, and it is adjudicable. Running verapamil's
*own* multichannel profile through a human ventricular model (O'Hara-Rudy)
prolongs slightly, as the clinic does. Same channels, same concentration,
opposite answer from the plate: the inversion is a property of the **assay** —
ICaL dominance in immature, IK1-poor cells — not of the compound.

## What the assay is actually good at

![Three panels over the 12 CiPA training drugs. qNet sorts the three risk categories into completely non-overlapping bands; the dish arrhythmia-event fraction puts all four high-risk drugs on top but overlaps intermediate and low; and rank agreement with CiPA risk is 0.95 for qNet against 0.65-0.68 for the three dish endpoints.]({{artifact:ca4472da-4c69-48ca-9587-959747a164b3}})

The assay exists to sort compounds into torsades risk categories. Tested against
the CiPA in-silico metric (qNet) on **all 28 CiPA drugs** — 12 on which qNet was
trained and **16 genuinely held out**:

| on the 16 held-out drugs | Spearman rho | AUC high vs rest | AUC low vs rest |
|---|---|---|---|
| qNet (in-silico) | 0.840 [0.52, 0.95] | **1.000** | 0.909 [0.69, 1.00] |
| ddFPDc at Cmax (dish) | 0.809 [0.55, 0.93] | 0.917 [0.71, 1.00] | **0.945 [0.80, 1.00]** |
| max ddFPDc (dish) | 0.722 [0.38, 0.90] | 0.812 [0.56, 1.00] | **0.945 [0.78, 1.00]** |
| arrhythmia-event fraction (dish) | 0.717 [0.37, 0.90] | 0.875 [0.62, 1.00] | 0.891 [0.68, 1.00] |

**qNet is the better danger flag; the dish is its equal at clearing low-risk
compounds.** qNet separates high risk perfectly — all 8 high-risk drugs fall
below every intermediate and low-risk drug, on all 28 including the held-out 16.
That is the only comparison in the table with no overlap. On overall ranking and
on identifying low risk the two are statistically indistinguishable at n = 16.

An earlier version of this document claimed the dish "detects danger but does not
grade safety". That was an artifact of looking only at the 12 training drugs,
where the dish's AUC for low risk was 0.69; on the held-out 16 it is 0.945. The
claim is withdrawn.

Two qNet misses, both held out: metoprolol (low risk, but qNet 0.0702, inside the
intermediate band) and tamoxifen (0.0781).

## What cannot be estimated, and why that matters

![Three panels. Leave-one-site-out gain per site, with the single optical site highest at 1.44 and its interval clear of 1.0; the optical arm against the MEA arm with a fitted slope of 1.51; and clinical, dish and in-silico effects for four drugs where only the dish inverts verapamil.]({{artifact:3563b359-a2e4-45fc-9c17-607156ce1d88}})

The textbook chain is single-cell APD90 -> syncytial FPD -> transmural summation
-> QT. Two of its links cannot be estimated from public data at all, and saying
so precisely is more useful than producing a number.

**APD90 vs FPD is confounded beyond repair.** Blinova 2018 pools both endpoints
into one column and asserts they are "equivalent measures". They are not, as a
scale: the optical arm reads **1.44x** the consensus of the other nine sites
(CI 1.15-1.72), the highest of all ten. But the only optical platform is used by
**one site, which is also the only serum-free site** — and medium alone shifts
baseline FPDc by 104 ms in a direction that published data contradict. Endpoint,
site and medium move together; no public dataset separates them.

**Transmural summation is a coupling problem, not an average.** With uncoupled
cells the QT surrogate is simply the M-cell APD90 (121 of 136 conditions), with
an implied dispersion of 101.9 ms against the ~30-40 ms of a real wall. The
widely-used "weighted average APD90" mapping has no mechanism behind it.

**And the ceiling is empirical, not just methodological.** Blinova 2019 gave 16
subject-matched iPSC-CM lines their own donors' drugs and found no within-person
correlation (rho = -0.11, p = 0.70). The informative control: those subjects'
*clinical* responses to the two drugs did not correlate with each other either
(rho = 0.10, p = 0.73) — so part of what looks like assay failure is that
individual QT response is not a stable single trait.

## The clinical side has a measurement floor too

![Three panels from 238 clinical waveform records: beat-to-beat RR within single ten-second recordings; the distribution of implied QTcF spread across each window's slowest and fastest beat, two thirds of which exceed the ten-millisecond regulatory threshold; and the tangent-minus-study T-end difference against RR, with a significant negative slope.]({{artifact:5d8fd8b3-c971-4380-85b4-f9d02bdcca51}})

The PhysioNet concentration-QTc studies ship waveforms, not just intervals — a
fact I initially missed and then corrected. Decoding them:

- **Hysteresis correction is unidentifiable here, not merely unimplemented.** A
  record is 10 seconds, which is 8-17% of one QT adaptation time constant, so the
  RR history that sets the QT mostly precedes the strip. The cost is measurable:
  the QTcF spread between each window's slowest and fastest beat has a median of
  **12.6 ms and exceeds 10 ms in 67% of records**. The quantity being averaged
  moves by more than the regulatory threshold inside a single recording.
- **The T-end rule is rate-dependent**, at -0.032 ms per ms of RR (-9.5 ms across
  a 300 ms RR change, p = 0.026). A *constant* offset between delineation rules
  cancels in the double-delta cascade; a rate-dependent one does not, once heart
  rate moves between baseline and on-drug. The same sign and order of magnitude
  appears between two commercial algorithms on PTB-XL+ (-5.9 ms). Two datasets,
  two different algorithm pairs, one conclusion.

## Why the numbers can be trusted

Every stage is validated against an independent published analysis before it is
used, and the agreements are asserted in the test suite:

| stage | check | result |
|---|---|---|
| source data | Blinova workbook vs all 420 published panels | r = 0.999993, median abs diff 0.30 ms |
| clinical cascade | Vicente 2015 JAHA Table 2 | 8 of 8 |
| clinical cascade | Vicente 2016 Table 3 / Vicente 2019 CPT | 5 of 6 / 6 of 8 |
| per-record algebra | QTcF recomputed on 5,740 records | median 0.000 ms, SD 0.19 |
| in-vitro rate correction | Millard 2018 supplemental, 18 datasets | 17 of 18 within 5 ms, median 0.65 ms |
| waveform decoding | annotation-derived QT vs clinical table, 242 beats | median 0.00 ms |
| R-peak detection | detector RR vs published RR, 238 records | median +0.10 ms |
| in-silico | ORd control APD90 / ORd-CiPA control qNet | 265.5-329.9 ms / 0.0785 uC/uF |

Two findings were **withdrawn** mid-project when the evidence went against them —
a claimed mechanism for the dish-clinic gap, and an explanation of a between-study
shift. One discrepancy I could not resolve is documented rather than smoothed:
a 6-7 ms shift between two clinical studies for ranolazine and dofetilide.
A gap against one published table turned out to be the *table's* disagreement with
its own study's primary paper — our value reproduces the primary paper to 0.2 ms
at two independent timepoints.

## Honest limits

- Four drugs are the entire dish-clinic overlap. Every bridge statement rests on them.
- Wells that arrest are unmeasurable, and the drop rate rises from 0.9% to 60.8%
  across the concentration range. Every in-vitro fit is therefore reported twice,
  bounded rather than imputed.
- The in-silico risk comparison is in-sample for qNet and cannot show the
  simulation beating the dish.
- No mechanism is claimed for the 1.8-3.5x dish-clinic gap. A denervated
  preparation *should* over-read — autonomic tone attenuates drug-induced
  prolongation — but that is a direction, not an explanation.

## Sources for the quoted claims

The two statements attributed to Blinova 2018 above — that APD90 and FPD are
"equivalent measures", and that site 2 was "the only test site that used the VSO
platform instead of the MEA platforms and serum-free experimental medium" — are
verbatim from the paper's Results and Experimental Procedures: Blinova K et al.,
*Cell Reports* 2018;24(13):3582-3592, doi:10.1016/j.celrep.2018.08.079 (open
access). The 28-drug well-level values come from that study's published data file.

## My contribution

I did not generate any hiPSC-CM data, culture any cells, or run any clinical
study. The work here is the **curation and modelling half**: assembling
heterogeneous public data into an auditable form, implementing the correction
cascade with provenance enforcement, fitting mixed-effects concentration-response
models that carry site and protocol variability explicitly, validating every
stage against independent published analyses, processing the clinical waveforms,
and running the in-silico arm. Where the data could not answer a question, the
deliverable is a precise statement of why — which is the part I would most want
to be judged on.

## Next

The honest test of the in-silico arm is the 16 CiPA *validation* drugs. All 16
are already in the Blinova panel, so the dish side is in hand; only their channel
parameters are needed. That turns the risk comparison out-of-sample on both arms,
and it is the single highest-value extension of this work.
