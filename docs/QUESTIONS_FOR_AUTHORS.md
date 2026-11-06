# Questions for the CiPA myocyte study authors

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

---

Source: Blinova K et al., *Cell Reports* 2018;24(13):3582-3592, doi:10.1016/j.celrep.2018.08.079.
Each question traces to a specific limit documented in `PROJECT_REPORT.md` sections 7 and 10.
