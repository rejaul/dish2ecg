# Does the T-end rule explain the between-study shift? No.

`dish2ecg` carries an unresolved 6-7 ms discrepancy: for ranolazine and
dofetilide, ECGRDVQ and ECGCIPA report different baseline- and placebo-corrected
QTcF at matched concentration, and shared-placebo effects account for only part
of it. A rate-dependent T-end rule was the leading candidate. It is not the cause.

## Method

Both studies ship median-beat waveforms with the investigators' own fiducial
annotations. 3,163 ECGRDVQ and 2,515 ECGCIPA median beats were decoded for the
ranolazine, dofetilide and placebo arms. Decoding is exact: annotation-derived QT
reproduces each study's clinical table to a **median of 0.00 ms** (ECGCIPA, 399
records checked; the same check passed earlier for ECGRDVQ).

An independent tangent (Lepeschkin) T-end was measured on VCGMAG with the same
code for both studies, and QTcF recomputed from it.

Because the baseline/placebo cascade is **linear in QT**, the change it would
produce is exactly the cascade applied to the per-record difference
`delta = QTcF_tangent - QTcF_study`. That avoids re-deriving the whole pipeline
and gives the answer without approximation.

## Result

Per record the two rules differ substantially, and differ *between* studies:
median **-4.5 ms** (ECGRDVQ) and **-7.1 ms** (ECGCIPA), SD 30 and 51 ms.

After baseline and placebo subtraction, almost all of that cancels:

| study | drug | movement in ddQTcF (ms) |
|---|---|---|
| ECGRDVQ | Ranolazine | -0.94 [-2.17, 1.05] |
| ECGRDVQ | Dofetilide | -1.83 [-4.71, -0.16] |
| ECGCIPA | Ranolazine | -0.61 [-4.92, 6.18] |
| ECGCIPA | Dofetilide | +0.25 [-3.37, 4.81] |

The **between-study** difference in that movement - the part that could explain
the shift - is **-0.32 ms for ranolazine and -2.08 ms for dofetilide**, against a
discrepancy of 6-7 ms. Delineation explains essentially none of it for ranolazine
and at most a third for dofetilide, with intervals that include zero.

## Why the residual is small, and what that corrects

The rate-dependence is real: the tangent-minus-study difference moves at
**-0.0317 ms per ms of RR**. Earlier documents quoted that as "-9.5 ms across a
300 ms RR change", which is a sensitivity, not an exposure. In these studies the
placebo-corrected within-subject RR change is only **-33 to +18 ms**, so the
predicted delineation effect is **+0.33 to +1.03 ms** - matching the observed
movement of -1.8 to +0.3 ms.

So the mechanism is confirmed but its magnitude in a thorough-QT design is about
1 ms, because drug and placebo arms have similar heart-rate trajectories and the
cascade subtracts what they share. The earlier statement that extraction method
is "a live contributor to the unexplained between-study shift" is **too strong
and is withdrawn**: it is a live contributor to *absolute* QT comparisons across
algorithms, and a ~1 ms contributor to a placebo-corrected one.

**The 6-7 ms shift remains unexplained.**
