# Papers carrying numbers or data for dish2ecg

Found by sweeping the FDA ECG group **by author and by dataset** rather than by drug or
biomarker. That change of strategy is what surfaces the analysis papers: 238 unique records
across nine queries, triaged down to the ten below.

## Used, and what each one settled

| year | author | journal | DOI | what it gives | access |
|---|---|---|---|---|---|
| 2015 | Vicente | J Am Heart Assoc 4:None | 10.1161/JAHA.114.001615 | Table 2: baseline- and placebo-corrected effect at Cmax for the four bridge drugs. **Gave the 8/8 ECGRDVQ validation.** | open (PMC4579946) |
| 2019 | Vicente | Clin Pharmacol Ther 105:943-953 | 10.1002/cpt.1303 | CiPA Phase I primary endpoints. Revealed the specification: ddJ-Tc at Cmax **on day 3**, ddQTc on day 1 for chloroquine. **Gave the 6/8 ECGCIPA validation.** | open (PMC6654598) |
| 2016 | Vicente | PLoS One 11:e0163619 | 10.1371/journal.pone.0163619 | Table 3: ECGDMMLD effects at Tmax. Gave the earlier 5/6 validation. | open (PMC5201270) |
| 2016 | Johannesen | PLoS One 11:e0166925 | 10.1371/journal.pone.0166925 | The automated J-Tpeak/Tpeak-Tend algorithm (FDA ecglib). Documents the measurement chain: 3x10 s ECGs per timepoint, median beats, vector-magnitude lead via the Guldenring transform. | open (PMC5201230) |
| 2017 | Couderc | J Electrocardiol 50:769-775 | 10.1016/j.jelectrocard.2017.08.025 | Multi-algorithm JTpc comparison. Confirms the form `JTpc = JTp/RR^a` and attributes it to Johannesen. PMC's text extraction drops the superscript, so the exponent still comes from my empirical recovery (0.5799). | open (PMC5699480) |

## Read: the within-person test, and its answer

| 2019 | Blinova | Clin Transl Sci 12:687-697 | 10.1111/cts.12674 | **READ. Answer: no within-person correlation.** 16 subject-matched iPSC-CM lines vs the same subjects' clinical ddQTc slopes: Spearman rho = -0.106 (p=0.70) dofetilide, -0.089 (p=0.75) moxifloxacin. The clinical responses to the two drugs did not correlate with each other either (rho=0.096, p=0.73). Bounds the package to drug-level use. | open (PMC6853144) |

## Wanted, and closed

These are the two primary reports of the studies the whole package rests on. Every
validation so far has been against a secondary analysis of them.

| 2015 | Johannesen | Clin Pharmacol Ther 99:214-23 | 10.1002/cpt.205 | ECGDMMLD primary. PMC5421403 exists but `get_full_text_article` returned 0 characters, so it is indexed but not retrievable. Most likely to explain the unresolved 6-7 ms gap on dofetilide-alone ddQTcF, because Vicente 2016 reports the value without the model specification. | open (PMC5421403) |
| 2014 | Johannesen | Clin Pharmacol Ther 96:549-558 | 10.1038/clpt.2014.155 | ECGRDVQ primary. An independent second check on the four bridge drugs. | closed; DOI did not resolve in PubMed's converter |

## Also surfaced, lower priority

- **Boulay 2024**, *J Pharmacol Toxicol Methods* - JTpc and TpTec as proarrhythmic markers.
- **Strauss 2017**, *Circulation* - genetic variant risk score against drug-induced QT in these studies.
- **Vicente 2015**, *J Electrocardiol* - sex differences in drug-induced repolarization, ECGRDVQ subgroups.
- **Darpo 2015**, *Clin Pharmacol Ther* / *Drug Saf* - the IQ-CSRC study, the design that replaced the thorough QT study.
- **Schocken 2018**, *J Pharmacol Toxicol Methods* (you supplied it) - media effects on hiPSC-CM readouts. Bears on in-vitro variability; not a data source.

## What is now validated

| corpus | source | agreement |
|---|---|---|
| ECGRDVQ | Vicente 2015 JAHA Table 2 | **8 of 8** |
| ECGCIPA | Vicente 2019 CPT, day-stratified | **6 of 8** |
| ECGDMMLD | Vicente 2016 PLoS ONE Table 3 | 5 of 6 |
| in-vitro rate correction | Millard 2018 Suppl. Table 1 | 17 of 18 within 5 ms |

Remaining disagreements, all reported rather than explained away: dofetilide-alone ddQTcF and
ddJ-Tpeakc in ECGDMMLD (6-7 ms low, specification not recoverable); verapamil ddQTcF in ECGCIPA
(5.7 ms low on day 1, 3.3 ms high on day 3 - a linear model on a curved concentration-response,
which the paper itself flags); and study 9 of the Millard table (22 ms, unexplained).