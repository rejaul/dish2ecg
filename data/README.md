# data/

This repository **does not redistribute third-party data**. It ships one derived
table of its own and tells you how to obtain the rest.

## Shipped here

| file | what it is |
|---|---|
| `cipa_concentration_lookup.csv` | Drug, concentration index 1-4, molar concentration, multiple of free Cmax, and Cmax, for all 28 drugs of the CiPA myocyte validation study. Read from the per-drug tables printed under the figure panels of the open-access supplement (Data S1-S3 of doi:10.1016/j.celrep.2018.08.079). Every row carries its source PDF and page, and a `conc_value_source` flag distinguishing a printed value from one reconstructed off the log ladder. |

A unit trap worth knowing about: those printed tables are headed "Drug
Concentration, M" but the values are **micromolar**. The micro sign is a
Symbol-font glyph that extracts as a private-use codepoint, so the unit reads as
plain "M". Dofetilide's printed Cmax of 0.002 is 2 nM, its known free clinical
Cmax; read as molar it would be 2 mM. The lookup stores mol/L.

## Fetched on demand

**Clinical ECG — PhysioNet.** `dish2ecg.io.load_physionet_clinical("ecgrdvq")`
and `("ecgdmmld")` download the clinical tables on first use, verify them against
the study's own `SHA256SUMS.txt`, and cache them under `~/.cache/dish2ecg`
(override with `DISH2ECG_CACHE`). Waveform directories are not downloaded.

- ECGRDVQ v1.0.0 — Johannesen L et al., *Clin Pharmacol Ther* 2014;96:549-558,
  doi:10.1038/clpt.2014.155
- ECGDMMLD v1.0.0 — Johannesen L et al., *Clin Pharmacol Ther* 2015;99:214-223,
  doi:10.1002/cpt.205
- ECGCIPA v1.0.0 — Vicente J et al., *Clin Pharmacol Ther* 2019;105:943-953,
  doi:10.1002/cpt.1303. Loaded with `dish2ecg.io.load_ecgcipa()`, which returns
  the four CDISC ADaM tables (`adsl`, `adeg`, `adpc`, `addm`) rather than one flat
  file. Watch the units: dofetilide is reported in pg/mL and every other analyte
  in ng/mL.

**In-vitro — the CiPA myocyte workbook.** `Blinova_etal_2018_data.xlsx` is the
well-level release of the FDA/HESI CiPA myocyte validation study, distributed
from the CiPA Data Resources page. Download it, then:

    python scripts/build_wells.py --workbook Blinova_etal_2018_data.xlsx

which joins it to the concentration lookup and writes
`data/cipa_myocyte_wells.csv`.

Cite: Blinova K et al., *Cell Rep* 2018;24:3582-3592,
doi:10.1016/j.celrep.2018.08.079.

**In-vitro validation reference — the CiPA MEA pilot.**
`scripts/validate_invitro_correction.py` downloads the supplementary archive of
Millard D et al., *Toxicol Sci* 2018;164:550-562, doi:10.1093/toxsci/kfy110
(PMC6061700) from Europe PMC and recovers Supplemental Table 1 — mean baseline
beat period, field potential duration and rate-corrected field potential duration
for 18 studies. Nothing is cached or redistributed.

## Two things to know before you model the in-vitro data

1. `ddFPDc` is **already** rate-, baseline- and vehicle-corrected by the source
   study. Do not run the `dish2ecg` cascade over it again.
2. Its missingness is **informative**. Of 2067 missing endpoints, 1887 are wells
   flagged with an arrhythmic event: a well in arrhythmia has no measurable field
   potential duration. Missingness rises from 0.9 % at the lowest concentration to
   60.8 % at the highest, so a complete-case slope is biased toward zero.

`drug_reference.csv` - one row per agent (15), with free Cmax, fraction unbound and a
per-cell source string; provenance and caveats in `docs/DRUG_REFERENCE_SOURCES.md`.
