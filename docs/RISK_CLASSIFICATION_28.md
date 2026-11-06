# The out-of-sample test: all 28 CiPA drugs, 16 held out

The earlier comparison (`RISK_CLASSIFICATION.md`) used the 12 CiPA *training*
drugs, on which the ORd-CiPA model and qNet were optimised, so qNet's performance
there was in-sample and no comparison against the dish was supportable.

The parameters for the remaining 16 turned out to be public after all, on the
FDA/CiPA repository's `Model-Validation-2018` branch: `newCiPA.csv` carries
7-channel IC50s, Hill coefficients, therapeutic concentrations and risk categories
for 29 drugs, and `hERG_fitting_results.zip` carries 2,000 bootstrap hERG binding
parameter sets for all 28. Those 28 are **exactly the Blinova panel**, so both
arms cover the same compounds.

## Handling the missing point estimates

The repository publishes an *optimal* hERG point fit only for the 12 training
drugs; for all 28 it publishes bootstrap samples. Marginal medians of those
samples are not a safe substitute — bepridil's median Kmax differs 5-fold from
its point fit — because Kmax and halfmax are strongly co-dependent.

So qNet was computed **per bootstrap sample** and summarised in output space:
median qNet as the point estimate, 2.5-97.5 percentiles as the interval,
20 samples per drug x 4 doses x 28 drugs = **2,240 simulations, 0 failures**.

That choice is validated where it can be: on the 12 training drugs, bootstrap-
median qNet matches the published point-fit qNet to within **0.0006 uC/uF for 11
of 12** (quinidine 0.0033), rank correlation 0.986.

## Result

qNet separates high risk **perfectly across all 28 drugs**: all 8 high-risk
compounds fall below every intermediate and low-risk compound (AUC 1.000),
including the 16 that were held out.

On the 16 held-out drugs:

| predictor | Spearman rho | AUC high vs rest | AUC low vs rest |
|---|---|---|---|
| qNet (in-silico) | 0.840 [0.52, 0.95] | **1.000** | 0.909 [0.69, 1.00] |
| ddFPDc at Cmax (dish) | 0.809 [0.55, 0.93] | 0.917 [0.71, 1.00] | **0.945 [0.80, 1.00]** |
| max ddFPDc (dish) | 0.722 [0.38, 0.90] | 0.812 [0.56, 1.00] | **0.945 [0.78, 1.00]** |
| arrhythmia-event fraction (dish) | 0.717 [0.37, 0.90] | 0.875 [0.62, 1.00] | 0.891 [0.68, 1.00] |

4,000-sample bootstrap intervals. **qNet is the better danger flag; the dish is
its equal at clearing low-risk compounds.** The perfect high-risk separation is
the only comparison in the table with no overlap; on ranking and on low-risk
identification the two are indistinguishable at this sample size.

Across all 28 drugs qNet leads on every metric (Spearman 0.899 [0.76, 0.94],
AUC high 1.000, AUC low 0.959 [0.85, 1.00] against a best dish value of 0.752,
0.919 and 0.912) — but 12 of those 28 are its training set, so the held-out
column above is the one to quote.

## What this corrects

The 12-drug analysis concluded the dish "detects danger but does not grade
safety", on an AUC for low risk of 0.69. On the held-out 16 that figure is 0.945.
**The claim was an artifact of the training subset and is withdrawn.**

## Where qNet fails

Two low-risk drugs sit inside the intermediate qNet band, both held out:
**metoprolol** (0.0702) and **tamoxifen** (0.0781). Neither is flagged as high
risk, so the regulatory consequence is a false alarm rather than a missed danger.
