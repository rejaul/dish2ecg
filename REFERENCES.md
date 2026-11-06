# References

Every record below was resolved against CrossRef on retrieval; author, year, journal and
pagination are as the registry returns them, not as recalled. Notes say what each paper is
used for in this package and whether its full text was obtained.

## Core — in vitro

- **Blinova K et al. (2018)** International Multisite Study of Human-Induced Pluripotent Stem Cell-Derived Cardiomyocytes for Drug Proarrhythmic Potential Assessment. *Cell Reports 24:3582-3592*. doi:10.1016/j.celrep.2018.08.079  
  Source of Blinova_etal_2018_data.xlsx; the 28-drug, 10-site hiPSC-CM study. ddFPDc already baseline-, vehicle- and Fridericia-corrected.
- **Millard D et al. (2018)** Cross-Site Reliability of Human Induced Pluripotent stem cell-derived Cardiomyocyte Based Safety Assays Using Microelectrode Arrays: Results from a Blinded CiPA Pilot Study. *Toxicological Sciences 164:550-562*. doi:10.1093/toxsci/kfy110  
  Cross-site reliability. Suppl. Table 1 is the 18-study target for the in-vitro rate-correction check (17/18 within 5 ms).
- **Blinova K et al. (2016)** Comprehensive Translational Assessment of Human-Induced Pluripotent Stem Cell Derived Cardiomyocytes for Evaluating Drug-Induced Arrhythmias. *Toxicological Sciences 155:234-247*. doi:10.1093/toxsci/kfw200  
  Earlier FDA hiPSC-CM translational study. Its clinical comparison table is categorical, not ms — a dead end for widening the bridge.
- **Blinova K et al. (2019)** Clinical Trial in a Dish: Personalized Stem Cell-Derived Cardiomyocyte Assay Compared With Clinical Trial Results for Two QT-Prolonging Drugs. *Clinical and Translational Science 12:687-697*. doi:10.1111/cts.12674  
  The subject-matched test. No within-person dish-to-clinic correlation; bounds the package to drug-level use.
- **Schocken D et al. (2018)** Comparative analysis of media effects on human induced pluripotent stem cell-derived cardiomyocytes in proarrhythmia risk assessment. *Journal of Pharmacological and Toxicological Methods 90:39-47*. doi:10.1016/j.vascn.2017.11.002  
  Media effects on hiPSC-CM proarrhythmia readouts. Attached, not yet mined.

## Core — clinical

- **Johannesen L et al. (2014)** Differentiating Drug-Induced Multichannel Block on the Electrocardiogram: Randomized Study of Dofetilide, Quinidine, Ranolazine, and Verapamil. *Clinical Pharmacology & Therapeutics 96:549-558*. doi:10.1038/clpt.2014.155  
  ECGRDVQ primary (dofetilide, quinidine, ranolazine, verapamil). CLOSED — still wanted.
- **Johannesen L et al. (2015)** Late sodium current block for drug-induced long QT syndrome: Results from a prospective clinical trial. *Clinical Pharmacology & Therapeutics 99:214-223*. doi:10.1002/cpt.205  
  ECGDMMLD primary (late-INa block). CLOSED — best candidate to explain the unresolved 6-7 ms dofetilide-alone gap.
- **Vicente J et al. (2019)** Assessment of Multi-Ion Channel Block in a Phase I Randomized Study Design: Results of the CiPA Phase I ECG Biomarker Validation Study. *Clinical Pharmacology & Therapeutics 105:943-953*. doi:10.1002/cpt.1303  
  ECGCIPA primary. Our day-stratified validation target: 6 of 8.
- **Vicente J et al. (2015)** Comprehensive T wave Morphology Assessment in a Randomized Clinical Study of Dofetilide, Quinidine, Ranolazine, and Verapamil. *Journal of the American Heart Association 4*. doi:10.1161/JAHA.114.001615  
  ECGRDVQ T-wave morphology; Table 2 is the 8/8 validation target. Note: its alpha=0.58 is for T-wave flatness, NOT for J-Tpeak.
- **Vicente J et al. (2016)** Electrocardiographic Biomarkers for Detection of Drug-Induced Late Sodium Current Block. *PLOS ONE 11:e0163619*. doi:10.1371/journal.pone.0163619  
  ECGDMMLD biomarkers; Table 3 is the 5/6 validation target.

## Methods — ECG intervals

- **Johannesen L et al. (2016)** Automated Algorithm for J-Tpeak and Tpeak-Tend Assessment of Drug-Induced Proarrhythmia Risk. *PLOS ONE 11:e0166925*. doi:10.1371/journal.pone.0166925  
  The automated J-Tpeak / Tpeak-Tend algorithm behind the PhysioNet interval columns.
- **Hnatkova K et al. (2019)** Heart Rate Correction of the J-to-Tpeak Interval. *Scientific Reports 9*. doi:10.1038/s41598-019-51491-4  
  Source of the LINEAR J-Tpeak correction JTpc = JTp + 0.150(1-RR). First author is Hnatkova, not Vicente.
- **Couderc J et al. (2017)** An evaluation of multiple algorithms for the measurement of the heart rate corrected JTpeak interval. *Journal of Electrocardiology 50:769-775*. doi:10.1016/j.jelectrocard.2017.08.025  
  Confirms the power-law form JTpc = JTp/RR^alpha; ECGCIPA's own JTPC fits alpha = 0.5799.

## Context — CiPA

- **Crumb WJ et al. (2016)** An evaluation of 30 clinical drugs against the comprehensive in vitro proarrhythmia assay (CiPA) proposed ion channel panel. *Journal of Pharmacological and Toxicological Methods 81:251-262*. doi:10.1016/j.vascn.2016.03.009  
  30-drug ion-channel panel. UNRETRIEVED — the tabulated fraction-unbound source for four CiPA drugs.
- **Li Z et al. (2017)** Improving the In Silico Assessment of Proarrhythmia Risk by Combining hERG (Human Ether-à-go-go-Related Gene) Channel-Drug Binding Kinetics and Multichannel Pharmacology. *Circulation: Arrhythmia and Electrophysiology 10*. doi:10.1161/CIRCEP.116.004628  
  In-silico proarrhythmia metric. UNRETRIEVED — second fraction-unbound source.
- **Colatsky T et al. (2016)** The Comprehensive in Vitro Proarrhythmia Assay (CiPA) initiative — Update on progress. *Journal of Pharmacological and Toxicological Methods 81:15-20*. doi:10.1016/j.vascn.2016.06.002  
  CiPA initiative update.
- **Strauss DG et al. (2019)** Comprehensive In Vitro Proarrhythmia Assay (CiPA) Update from a Cardiac Safety Research Consortium / Health and Environmental Sciences Institute / FDA Meeting. *Therapeutic Innovation & Regulatory Science 53:519-525*. doi:10.1177/2168479018795117  
  CiPA update from the CSRC/HESI/FDA meeting.

## Added with the research report

Every record resolved against CrossRef on retrieval; print years follow the journal volume where CrossRef records an online-first year.

### Clinical — J-Tpeak

- Johannesen L et al. (2014) Improving the Assessment of Heart Toxicity for All New Drugs Through Translational Regulatory Science. *Clinical Pharmacology & Therapeutics* 95:501-508. doi:10.1038/clpt.2013.238  
  Attributed by Vicente 2019 as the origin of the J-Tpeakc exponent 0.58. NOT READ here.
- Hnatkova K et al. (2017) Heart rate dependency of JT interval sections. *Journal of Electrocardiology* 50:814-824. doi:10.1016/j.jelectrocard.2017.08.005  
  Rate dependence of JT interval sections; related to the J-Tpeak correction. Not mined.

### Clinical — design

- Vicente J et al. (2018) Mechanistic Model-Informed Proarrhythmic Risk Assessment of Drugs: Review of the “CiPA” Initiative and Design of a Prospective Clinical Validation Study. *Clinical Pharmacology & Therapeutics* 103:54-66. doi:10.1002/cpt.896  
  CiPA review and design of the prospective ECG validation study that became ECGCIPA.

### Cited in the research report

- Malik M et al. (2019) Implications of Individual QT/RR Profiles—Part 1: Inaccuracies and Problems of Population-Specific QT/Heart Rate Corrections. *Drug Safety* 42:401-414. doi:10.1007/s40264-018-0736-1  
  Cited in docs/RESEARCH_REPORT.md.

### Context — CiPA and in-vitro platforms

- Li Z et al. (2019) Assessment of an In Silico Mechanistic Model for Proarrhythmia Risk Prediction Under the CiPA Initiative. *Clinical Pharmacology & Therapeutics* 105:466-475. doi:10.1002/cpt.1184  
  In-silico mechanistic model for proarrhythmia risk under CiPA.
- Dutta S et al. (2017) Optimization of an In silico Cardiac Cell Model for Proarrhythmia Risk Assessment. *Frontiers in Physiology* 8. doi:10.3389/fphys.2017.00616  
  Optimisation of the in-silico cardiac cell model; source of the free-Cmax convention for the CiPA training drugs.
- Sager PT et al. (2014) Rechanneling the cardiac proarrhythmia safety paradigm: A meeting report from the Cardiac Safety Research Consortium. *American Heart Journal* 167:292-300. doi:10.1016/j.ahj.2013.11.004  
  Meeting report that launched CiPA.
- Fermini B et al. (2016) A New Perspective in the Field of Cardiac Safety Testing through the Comprehensive In Vitro Proarrhythmia Assay Paradigm. *SLAS Discovery* 21:1-11. doi:10.1177/1087057115594589  
  CiPA paradigm overview.
- Kilfoil P et al. (2021) Characterization of a high throughput human stem cell cardiomyocyte assay to predict drug-induced changes in clinical electrocardiogram parameters. *European Journal of Pharmacology* 912:174584. doi:10.1016/j.ejphar.2021.174584  
  High-throughput hiPSC-CM assay for torsades prediction.
- Ernault AC et al. (2024) Interpretation of field and LEAP potentials recorded from cardiomyocyte monolayers. *American Journal of Physiology-Heart and Circulatory Physiology* 326:H800-H811. doi:10.1152/ajpheart.00463.2023  
  Interpretation of field and LEAP potentials from monolayers; the physics behind the FPD waveform.
- Guerrelli D et al. (2024) hiPSC-CM electrophysiology: impact of temporal changes and study parameters on experimental reproducibility. *American Journal of Physiology-Heart and Circulatory Physiology* 327:H12-H27. doi:10.1152/ajpheart.00631.2023  
  hiPSC-CM MEA reproducibility; well exclusion criteria and temporal drift.

### Clinical — concentration-QTc

- Garnett C et al. (2016) Operational Characteristics of Linear Concentration-QT Models for Assessing QTc Interval in the Thorough QT and Phase I Clinical Studies. *Clinical Pharmacology & Therapeutics* 100:170-178. doi:10.1002/cpt.361  
  Operational characteristics of linear concentration-QT models.
- Garnett C et al. (2018) Scientific white paper on concentration-QTc modeling. *Journal of Pharmacokinetics and Pharmacodynamics* 45:383-397. doi:10.1007/s10928-017-9558-5  
  Scientific white paper on concentration-QTc modelling; the linear mixed-effects model this package implements.

### Subject-level studies

- Stillitano F et al. (2017) Modeling susceptibility to drug-induced long QT with a panel of subject-specific induced pluripotent stem cells. *eLife* 6. doi:10.7554/elife.19406  
  Subject-specific iPSC-CMs from extreme sotalol responders; in-vitro responses correlated with in-vivo. Prior positive finding.
- Shinozawa T et al. (2017) Recapitulation of Clinical Individual Susceptibility to Drug-Induced QT Prolongation in Healthy Subjects Using iPSC-Derived Cardiomyocytes. *Stem Cell Reports* 8:226-234. doi:10.1016/j.stemcr.2016.12.014  
  Ten healthy volunteers, moxifloxacin; QT correlated with FPD at clinically relevant concentrations. Prior positive finding.

### Open software (screened for novelty)

- Dunham CS et al. (2022) Cardio PyMEA: A user-friendly, open-source Python application for cardiomyocyte microelectrode array analysis. *PLOS ONE* 17:e0266647. doi:10.1371/journal.pone.0266647  
  Cardio PyMEA, open-source MEA analysis. Signal processing only.
- Lee J et al. (2024) CardioMEA: comprehensive data analysis platform for studying cardiac diseases and drug responses. *Frontiers in Physiology* 15. doi:10.3389/fphys.2024.1472126  
  CardioMEA, open-source HD-MEA analysis. Signal processing only.

### Data and software

- Goldberger AL et al. (2000) PhysioBank, PhysioToolkit, and PhysioNet. *Circulation* 101. doi:10.1161/01.cir.101.23.e215  
  PhysioNet.
- Seabold S et al. (2010) Statsmodels: Econometric and Statistical Modeling with Python. *Proceedings of the Python in Science Conference*:92-96. doi:10.25080/majora-92bf1922-011  
  statsmodels.
- Virtanen P et al. (2020) SciPy 1.0: fundamental algorithms for scientific computing in Python. *Nature Methods* 17:261-272. doi:10.1038/s41592-019-0686-2  
  SciPy.
- Harris CR et al. (2020) Array programming with NumPy. *Nature* 585:357-362. doi:10.1038/s41586-020-2649-2  
  NumPy.
- McKinney W (2010) Data Structures for Statistical Computing in Python. *Proceedings of the Python in Science Conference*:56-61. doi:10.25080/majora-92bf1922-00a  
  pandas.

### Statistics

- Little R et al. (2019) Statistical Analysis with Missing Data, Third Edition. *Wiley Series in Probability and Statistics*. doi:10.1002/9781119482260  
  Missing-data theory: missing completely at random, at random, not at random.
- Pinheiro JC et al. (2000) Mixed-Effects Models in S and S-PLUS. *Statistics and Computing*. doi:10.1007/978-1-4419-0318-1  
  Mixed-effects models.
- Kenward MG et al. (1997) Small Sample Inference for Fixed Effects from Restricted Maximum Likelihood. *Biometrics* 53:983. doi:10.2307/2533558  
  Kenward-Roger degrees of freedom; the df rule common in statistical analysis plans. This package uses t on n_groups-1 instead.

## Public ECG corpora assessed (appendix 8.1 of the research report)

Diagnosis and arrhythmia corpora; none carries drug exposure. PTB-XL+ supplied the
interval-extraction measurement in uncertainty 9b.

- Wagner P et al. (2020) PTB-XL, a large publicly available electrocardiography dataset. *Scientific Data* 7. doi:10.1038/s41597-020-0495-6
- Strodthoff N et al. (2023) PTB-XL+, a comprehensive electrocardiographic feature dataset. *Scientific Data* 10. doi:10.1038/s41597-023-02153-8
- Moody GB et al. (2001) The impact of the MIT-BIH Arrhythmia Database. *IEEE Engineering in Medicine and Biology Magazine* 20:45-50. doi:10.1109/51.932724
- Clifford G et al. (2017) AF Classification from a Short Single Lead ECG Recording: the Physionet Computing in Cardiology Challenge 2017. *Computing in Cardiology Conference (CinC)* 44. doi:10.22489/cinc.2017.065-469
- Perez Alday EA et al. (2020) Classification of 12-lead ECGs: the PhysioNet/Computing in Cardiology Challenge 2020. *Physiological Measurement* 41:124003. doi:10.1088/1361-6579/abc960
- Reyna MA et al. (2021) Will Two Do? Varying Dimensions in Electrocardiography: The PhysioNet/Computing in Cardiology Challenge 2021. *2021 Computing in Cardiology (CinC)*:1-4. doi:10.23919/cinc53138.2021.9662687
- Liu F et al. (2018) An Open Access Database for Evaluating the Algorithms of Electrocardiogram Rhythm and Morphology Abnormality Detection. *Journal of Medical Imaging and Health Informatics* 8:1368-1373. doi:10.1166/jmihi.2018.2442
- Bousseljot R et al. (1995) Nutzung der EKG-Signaldatenbank CARDIODAT der PTB über das Internet. *Biomedizinische Technik/Biomedical Engineering*:317-318. doi:10.1515/bmte.1995.40.s1.317
- Zhang Jw et al. (2010) Chinese Cardiovascular Disease Database (CCDD) and Its Management Tool. *2010 IEEE International Conference on BioInformatics and BioEngineering*:66-72. doi:10.1109/bibe.2010.19
- Pilia N et al. (2021) ECGdeli - An open source ECG delineation toolbox for MATLAB. *SoftwareX* 13:100639. doi:10.1016/j.softx.2020.100639

## Translation-review literature (docs/TRANSLATION_REVIEW.md)

Retrieved in a November 2026 sweep of 28 PubMed queries plus targeted lookups; each resolved
against CrossRef. Not a systematic review.

- Antzelevitch C et al. (2007) Does Tpeak-Tend provide an index of transmural dispersion of repolarization?. *Heart Rhythm* 4:1114-1116. doi:10.1016/j.hrthm.2007.05.028
- Braam SR et al. (2013) Repolarization reserve determines drug responses in human pluripotent stem cell derived cardiomyocytes. *Stem Cell Research* 10:48-56. doi:10.1016/j.scr.2012.08.007
- Chang KC et al. (2017) Uncertainty Quantification Reveals the Importance of Data Variability and Experimental Design Considerations for in Silico Proarrhythmia Risk Assessment. *Frontiers in Physiology* 8. doi:10.3389/fphys.2017.00917
- Diedrich A et al. (2002) Modulation of QT Interval During Autonomic Nervous System Blockade in Humans. *Circulation* 106:2238-2243. doi:10.1161/01.cir.0000035241.76918.6c
- Gaur N et al. (2020) Validation of quantitative measure of repolarization reserve as a novel marker of drug induced proarrhythmia. *Journal of Molecular and Cellular Cardiology* 145:122-132. doi:10.1016/j.yjmcc.2020.04.019
- Goversen B et al. (2018) The immature electrophysiological phenotype of iPSC-CMs still hampers in vitro drug screening: Special focus on I K1. *Pharmacology & Therapeutics* 183:127-136. doi:10.1016/j.pharmthera.2017.10.001
- Hayes HB et al. (2019) Novel method for action potential measurements from intact cardiac monolayers with multiwell microelectrode array technology. *Scientific Reports* 9. doi:10.1038/s41598-019-48174-5
- Hazewinkel A et al. (2022) Sensitivity to missing not at random dropout in clinical trials: Use and interpretation of the trimmed means estimator. *Statistics in Medicine* 41:1462-1481. doi:10.1002/sim.9299
- Huo J et al. (2018) Sex-Related Differences in Drug-Induced QT Prolongation and Torsades de Pointes: A New Model System with Human iPSC-CMs. *Toxicological Sciences*. doi:10.1093/toxsci/kfy239
- Kanda Y et al. (2018) Development of torsadogenic risk assessment using human induced pluripotent stem cell-derived cardiomyocytes: Japan iPS Cardiac Safety Assessment (JiCSA) update. *Journal of Pharmacological Sciences* 138:233-239. doi:10.1016/j.jphs.2018.10.010
- Kettlewell S et al. (2019) Dynamic clamping human and rabbit atrial calcium current: narrowing ICaL window abolishes early afterdepolarizations. *The Journal of Physiology* 597:3619-3638. doi:10.1113/JP277827
- Kisitu J (2020) Chemical concentrations in cell culture compartments (C5) - free concentrations. *ALTEX*. doi:10.14573/altex.2008251
- Leishman DJ (2020) Improving prediction of torsadogenic risk in the CiPA in silico model by appropriately accounting for clinical exposure. *Journal of Pharmacological and Toxicological Methods* 101:106654. doi:10.1016/j.vascn.2019.106654
- Li M et al. (2017) Overexpression of KCNJ2 in induced pluripotent stem cell-derived cardiomyocytes for the assessment of QT-prolonging drugs. *Journal of Pharmacological Sciences* 134:75-85. doi:10.1016/j.jphs.2017.05.004
- Patel D et al. (2019) Assessment of Proarrhythmic Potential of Drugs in Optogenetically Paced Induced Pluripotent Stem Cell-Derived Cardiomyocytes. *Toxicological Sciences* 170:167-179. doi:10.1093/toxsci/kfz076
- Pfeiffer ER et al. (2016) Specific prediction of clinical QT prolongation by kinetic image cytometry in human stem cell derived cardiomyocytes. *Journal of Pharmacological and Toxicological Methods* 81:263-273. doi:10.1016/j.vascn.2016.04.007
- Qu Y et al. (2025) Delayed Repolarization Caused by                     hERG                     Block With Different Drug Modalities Can Be Detected in Stem Cell-Derived Cardiomyocytes: Incubation Time Matters. *Clinical and Translational Science* 18. doi:10.1111/cts.70283
- Rast G et al. (2016) Influence of field potential duration on spontaneous beating rate of human induced pluripotent stem cell-derived cardiomyocytes: Implications for data analysis and test system selection. *Journal of Pharmacological and Toxicological Methods* 82:74-82. doi:10.1016/j.vascn.2016.08.002
- SMITH AH et al. (2007) Autonomic Tone Attenuates Drug-Induced QT Prolongation. *Journal of Cardiovascular Electrophysiology* 18:960-964. doi:10.1111/j.1540-8167.2007.00901.x
- Schmid C et al. (2020) Characterization of iCell cardiomyocytes using single-cell RNA-sequencing methods. *Journal of Pharmacological and Toxicological Methods* 106:106915. doi:10.1016/j.vascn.2020.106915
- Seibertz F et al. (2023) Electrophysiological and calcium-handling development during long-term culture of human-induced pluripotent stem cell-derived cardiomyocytes. *Basic Research in Cardiology* 118. doi:10.1007/s00395-022-00973-0
- Soepriatna AH et al. (2023) Action potential metrics and automated data analysis pipeline for cardiotoxicity testing using optically mapped hiPSC-derived 3D cardiac microtissues. *PLOS ONE* 18:e0280406. doi:10.1371/journal.pone.0280406
- TeBay C et al. (2022) Metabolic and electrolyte abnormalities as risk factors in drug-induced long QT syndrome. *Biophysical Reviews* 14:353-367. doi:10.1007/s12551-022-00929-7
- Vargas HM et al. (2020) Time for a Fully Integrated Nonclinical-Clinical Risk Assessment to Streamline QT Prolongation Liability Determinations: A Pharma Industry Perspective. *Clinical Pharmacology & Therapeutics* 109:310-318. doi:10.1002/cpt.2029
- Vicente J et al. (2015) Sex differences in drug-induced changes in ventricular repolarization. *Journal of Electrocardiology* 48:1081-1087. doi:10.1016/j.jelectrocard.2015.08.004
- Wei F et al. (2020) Effects of Electrical Stimulation on hiPSC-CM Responses to Classic Ion Channel Blockers. *Toxicological Sciences* 174:254-265. doi:10.1093/toxsci/kfaa010
- Wei F et al. (2022) Effects of Serum and Compound Preparation Methods on Delayed Repolarization Evaluation With Human iPSC-CMs. *Toxicological Sciences* 188:48-61. doi:10.1093/toxsci/kfac043
- Wu WW et al. (2023) ICH S7B In Vitro Assays Do Not Address Mechanisms of QTC Prolongation for Peptides and Proteins - Data in Support of Not Needing Dedicated QTC Studies. *Clinical Pharmacology & Therapeutics* 114:1332-1341. doi:10.1002/cpt.3047
- Yan GX et al. (1998) Characteristics and Distribution of M Cells in Arterially Perfused Canine Left Ventricular Wedge Preparations. *Circulation* 98:1921-1927. doi:10.1161/01.cir.98.18.1921
- da Rocha AM et al. (2019) Detection of Drug-Induced Torsades de Pointes Arrhythmia Mechanisms Using hiPSC-CM Syncytial Monolayers in a High-Throughput Screening Voltage Sensitive Dye Assay. *Toxicological Sciences* 173:402-415. doi:10.1093/toxsci/kfz235
