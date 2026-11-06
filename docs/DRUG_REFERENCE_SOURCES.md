# `drug_reference.csv` — provenance and coverage

Built 2026-11-02. One row per agent that appears in ECGRDVQ (SCR-002) or ECGDMMLD (SCR-003),
plus the 12 CiPA training drugs — 15 rows after taking the union (dofetilide, quinidine,
ranolazine, verapamil, mexiletine and diltiazem belong to both sets). `MOXI.M2` is carried as
a row because it is a measured analyte in ECGDMMLD, not because it is a dosed agent.

Every cell carries its own source: the `source_for_each_value` column holds a JSON object keyed
by column name. A value that could not be sourced is **empty**, with source `"not found"` —
nothing in this table is inferred from a drug class, a structural analogue, or a guess.

## Columns

| column | unit / coding | how it was obtained |
|---|---|---|
| `drug_name` | lowercase free base | normalised name; the ECGRDVQ `EXTRT` strings carry salts (`Quinidine Sulph`, `Verapamil HCL`) and the ECGDMMLD columns are three/four-letter codes — both join to this key by hand in the loader |
| `inchikey` | standard InChIKey | PubChem PUG-REST through the `chemistry` connector, name lookup, top hit; free base / parent stereoisomer |
| `unii` | FDA UNII, 10 chars | PubChem synonym of the form `UNII-<code>` for the same CID |
| `pubchem_cid` | integer | PubChem name query; the top hit of `n_cids_total` matches |
| `chembl_id` | ChEMBL molecule id | `chembl` connector `compound_search`, accepting only an exact `pref_name` match to the free base (salt entries such as `VERAPAMIL HYDROCHLORIDE` were rejected in favour of the parent) |
| `total_cmax_molar` | mol/L | **derived**, `free_cmax_molar / fraction_unbound` — not an independently measured value |
| `free_cmax_molar` | mol/L | `AP_simulation/data/CiPA_training_drugs.csv` column `therapeutic` (nM), FDA/CiPA commit `539ff80d6616895de5d943e7bfb77b346d195db4` |
| `fraction_unbound` | fraction, 0–1 | `1 −` plasma-protein-binding percentage quoted verbatim from an FDA SPL label retrieved through openFDA; where the label gives a range the midpoint is used and the range is quoted in the source string |
| `fu_low`, `fu_high`, `fu_basis` | fraction; text | added 2026-11-04; populated for the four drugs above only. Other rows carry their range inside `source_for_each_value` |
| `cipa_risk_category` | high / intermediate / low | `CiPA_training_drugs.csv` column `CiPA` |

## Why `free_cmax_molar` is the *free* concentration

The FDA/CiPA repository labels the column only as "therapeutic concentration (nM)"
(`AP_simulation/README.md`). Two open-access CiPA papers state the quantity explicitly:
Dutta et al. 2017, *Front Physiol* 8:616 defines Cmax as "free maximum plasma clinical drug
exposures", and Li et al. 2019, *Clin Pharmacol Ther* 105:466 refers to the CiPA concentration
range as multiples of the "maximum free plasma concentration (Cmax)". Both were retrieved in
full text for this check. The `therapeutic` column is therefore read as free Cmax.

## Risk-category coding

`CiPA_training_drugs.csv` encodes risk as 2 / 1 / 0. The mapping 2 = high, 1 = intermediate,
0 = low was not taken on faith: it was confirmed against the risk headings of Supplemental
Table 2 of Blinova et al. 2018, *Cell Rep* 24:3582 (PMC6226030), which groups the same drugs
under "High TdP risk" (bepridil, d,l-sotalol, dofetilide, quinidine), "Intermediate TdP risk"
(chlorpromazine, cisapride, ondansetron, terfenadine) and "Low TdP risk" (diltiazem,
mexiletine, ranolazine, verapamil). The two sources agree on all 12 drugs.

## Cross-check of free Cmax against Blinova 2018 Table S2

Blinova's Table S2 lists its own Cmax (uM) for the same 12 drugs. Eight agree with the
CiPA training file to within ~10 % (dofetilide 0.002 vs 0.002, ranolazine 1.948 vs 1.9482,
bepridil 0.032 vs 0.033, sotalol 15 vs 14.69, quinidine 3 vs 3.237, diltiazem 0.128 vs 0.122,
chlorpromazine 0.0345 vs 0.038, cisapride 0.00258 vs 0.0026). Four differ materially (Blinova
value first): verapamil 0.045 vs 0.081, mexiletine 2.5 vs 4.129, ondansetron 0.372 vs 0.139,
terfenadine 0.000286 vs 0.004 uM. **The table carries the CiPA training-file value**, because
that is the value the published Hill/IC50 and qNet machinery is calibrated against; the Blinova
value is recorded here only as a cross-check. A dish-to-clinic comparison that places the
Blinova well-level data on a fold-Cmax axis must use Blinova's Cmax, not this column, for those
four drugs.

## Coverage, honestly

| column | populated | missing |
|---|---|---|
| `inchikey` | 15/15 | — |
| `pubchem_cid` | 15/15 | — |
| `unii` | 14/15 | bepridil |
| `chembl_id` | 14/15 | moxifloxacin_N_sulfate_M2 |
| `free_cmax_molar` | 12/15 | lidocaine, moxifloxacin, moxifloxacin_N_sulfate_M2 |
| `cipa_risk_category` | 12/15 | lidocaine, moxifloxacin, moxifloxacin_N_sulfate_M2 |
| `fraction_unbound` | 14/15 | moxifloxacin_N_sulfate_M2 |
| `total_cmax_molar` | 12/15 | lidocaine, moxifloxacin, moxifloxacin_N_sulfate_M2 |

Specifics:

1. **bepridil has no UNII here.** PubChem CID 2351 carries 145 synonyms, none in `UNII-<code>`
   form; ChEBI:3061 has CAS and registry cross-references but no FDA SRS link. One bare
   ten-character token (`755BO701MA`) appears in the synonym list and is plausibly the UNII,
   but nothing in the retrieved records labels it as such, so the cell is left empty rather
   than filled with an unlabelled string.
2. **The moxifloxacin M2 metabolite has no ChEMBL id.** `compound_search(name="moxifloxacin
   N-sulfate")` returns zero compounds. PubChem resolves it (CID 71750858, UNII 9EI1SIA8XV).
3. **lidocaine, moxifloxacin and the M2 metabolite have no free Cmax and no CiPA risk class.**
   They are not CiPA training drugs, so they are absent from `CiPA_training_drugs.csv`, and
   they are absent from Blinova's 28-drug Table S2. Li et al. 2019 (the CiPA validation-drug
   paper) was retrieved in full text, but its drug-level table is carried in figures that the
   text extraction does not contain, so no risk class could be read from it. Searching for
   these values in a non-CiPA source would mix concentration conventions, so the cells are empty.
4. **Fraction unbound for chlorpromazine, bepridil, cisapride and terfenadine (added 2026-11-04).**
   Earlier versions left these four empty and stated that Crumb et al. 2016 and Li et al. 2017
   tabulate protein binding for this drug set. **That was wrong**: both PDFs were later supplied and
   searched, and neither contains a protein-binding or total-Cmax value for these drugs (Crumb 2016
   Table 1 lists only free Cmax and the four test concentrations; the Li 2017 hits are about hERG
   binding kinetics, not plasma binding). The values now come from product labels and measured
   literature instead, each stored with its source string, a range, and a `fu_basis`:

   | drug | fu (point) | range | basis |
   |---|---|---|---|
   | bepridil | 0.005 | 0.003-0.010 | label >99% bound; measured >99.7% (Pritchard 1985); saturable AAG binding |
   | chlorpromazine | 0.055 | 0.03-0.08 | 92-97% bound (label-derived monograph); compendium 96% |
   | cisapride | 0.0225 | 0.020-0.025 | label 97.5-98% bound |
   | terfenadine | 0.03 | single value | Seldane label 97% bound (1995) |

   Caveats that matter for use: (a) none of the four has a current FDA SPL, so the label text comes
   from reproductions on drug-information sites, not from openFDA; chlorpromazine is not a label
   value at all but an AHFS-derived monograph. (b) Bepridil's binding is saturable, so its fu is
   concentration-dependent and 0.005 is a bounded estimate (geometric mean of the two bounds), not a
   measurement; it agrees with the fu implied by Crumb free Cmax over the label's total Cmax (0.0050).
   (c) For cisapride the published peak totals (150-300 nM) imply fu 0.009-0.017, below the label
   range - a mild tension, recorded not resolved. (d) fu is not a constant: basic drugs bind
   alpha-1-acid glycoprotein, which rises in inflammation and malignancy.

   **Terfenadine's free Cmax has a different meaning from the others.** The table's CiPA value
   (4 nM) is 13x Crumb's and Blinova's (0.3 nM). With fu 0.03 the two imply total concentrations of
   ~63 and ~4.7 ng/mL; the reported normal-dose mean peak is ~1.5 ng/mL. A Redfern et al. 2003
   (Cardiovasc Res 58:32-45) excerpt states the 4 nM figure is a concentration attained in the
   presence of a P450 inhibitor - read from a search-result excerpt of the OUP PDF; the full text is
   closed and was not retrieved. So the 13x divergence is consistent with the CiPA value being an
   interaction-condition exposure, and a fold-Cmax axis for terfenadine means something different
   depending on which Cmax is used.

   None of these four drugs is in the dish-clinic bridge (no clinical concentration-QTc data), so
   adding their fu changes no published result of this repository.

5. **`total_cmax_molar` is a derived quantity, not a measurement.** It is `free / fu`, so it
   exists only for the 8 drugs that have both, and it inherits the imprecision of a label
   protein-binding range quoted as a midpoint. For quinidine (fu 0.16 from an 80–88 % binding
   range) the derived total, 20.2 µM, sits above the usual therapeutic total plasma range;
   treat derived totals as order-of-magnitude anchors, not as reference values.

## Citations for the values in this table

* FDA/CiPA repository, `AP_simulation/data/CiPA_training_drugs.csv`, commit `539ff80d6616895de5d943e7bfb77b346d195db4`
  (https://github.com/FDA/CiPA) — free Cmax and TdP risk category.
* Dutta, S., Chang, K.C., Beattie, K.A., et al. (2017). Optimization of an in silico cardiac
  cell model for proarrhythmia risk assessment. *Front Physiol* 8:616.
  doi:10.3389/fphys.2017.00616 — definition of Cmax as free plasma concentration.
* Li, Z., Ridder, B.J., Han, X., et al. (2019). Assessment of an in silico mechanistic model
  for proarrhythmia risk prediction under the CiPA initiative. *Clin Pharmacol Ther* 105:466.
  doi:10.1002/cpt.1184 — same definition, validation-drug context.
* Blinova, K., Dang, Q., Millard, D., et al. (2018). International multisite study of
  human-induced pluripotent stem cell-derived cardiomyocytes for drug proarrhythmic potential
  assessment. *Cell Rep* 24:3582. doi:10.1016/j.celrep.2018.08.079, PMC6226030,
  Supplemental Table 2 — independent Cmax column and risk-category headings.
* FDA Structured Product Labels via openFDA (`drug-regulatory` connector, `search_drug_labels`,
  `clinical_pharmacology` section) — plasma protein binding; the quoted sentence for each drug
  is stored in `source_for_each_value`.
* PubChem PUG-REST and ChEMBL via the `chemistry` and `chembl` connectors — CID, InChIKey,
  UNII synonym, ChEMBL id.
