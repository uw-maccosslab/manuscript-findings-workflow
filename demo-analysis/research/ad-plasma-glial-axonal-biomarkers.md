---
# Research-finding document — external knowledge, structured like a finding but for the literature/tools.
# Spec: doc 04.4. Stored in research/<slug>.md. A research finding without verified references is NOT accepted.

topic: "Glial (astrocytic), axonal, and synaptic neurodegeneration protein biomarkers of Alzheimer's disease measurable in blood plasma (NEFL, GFAP, NRGN, CHI3L1/YKL-40; GAP43 and SNAP25 as CSF-only comparators)"
type: research-publications
slug: ad-plasma-glial-axonal-biomarkers
created: 2026-09-15
updated: 2026-09-15

status: reviewed
reviewed_by: research-reviewer (claude@ogdb.com)
reviewed_date: 2026-09-15

# Normalized domain-entity references this research is about (canonical IDs = UniProt accessions, verified
# directly against UniProt; label = HGNC gene symbol).
entities:
  - { db: uniprot, id: "P07196", label: "NEFL" }
  - { db: uniprot, id: "P14136", label: "GFAP" }
  - { db: uniprot, id: "Q92686", label: "NRGN" }
  - { db: uniprot, id: "P36222", label: "CHI3L1" }
  - { db: uniprot, id: "P17677", label: "GAP43" }
  - { db: uniprot, id: "P60880", label: "SNAP25" }

# (required, non-empty) Every external claim traces to a reference here. Each is fact-checked by the
# research-reviewer: it must EXIST and SUPPORT the claim attributed to it.
references:
  - id: "doi:10.1001/jamaneurol.2016.6117"
    type: doi
    doi: "10.1001/jamaneurol.2016.6117"
    pmid: "28346578"
    citation: "Mattsson N, Andreasson U, Zetterberg H, Blennow K; Alzheimer's Disease Neuroimaging Initiative. Association of Plasma Neurofilament Light With Neurodegeneration in Patients With Alzheimer Disease. JAMA Neurol. 2017;74(5):557-566."
    claim: "In the ADNI cohort (193 controls, 197 MCI, 180 AD dementia), plasma NfL (measured by an ultrasensitive assay) was higher in MCI (mean 42.8 ng/L) and AD dementia (mean 51.0 ng/L) than controls (mean 34.7 ng/L, P<.001), discriminated AD dementia from controls with AUC 0.87 (comparable to CSF biomarkers), correlated with CSF NfL (Spearman rho=0.59), and was associated with cognitive decline, AD-related atrophy, and brain hypometabolism."
    verified: true
    verified_by: research-reviewer (claude@ogdb.com)

  - id: "doi:10.1038/s41591-018-0304-3"
    type: doi
    doi: "10.1038/s41591-018-0304-3"
    pmid: "30664784"
    citation: "Preische O, Schultz SA, Apel A, et al; Dominantly Inherited Alzheimer Network. Serum neurofilament dynamics predicts neurodegeneration and clinical progression in presymptomatic Alzheimer's disease. Nat Med. 2019;25(2):277-283."
    claim: "In the Dominantly Inherited Alzheimer Network (autosomal-dominant AD), serum NfL (n=405, ultrasensitive immunoassay) correlated with CSF NfL, was elevated at presymptomatic stages, and its longitudinal rate of change discriminated mutation carriers from non-carriers ~16.2 years before estimated symptom onset (vs 6.8 years for cross-sectional NfL), tracking cortical thinning."
    verified: true
    verified_by: research-reviewer (claude@ogdb.com)

  - id: "doi:10.1515/cclm-2015-1195"
    type: doi
    doi: "10.1515/cclm-2015-1195"
    pmid: "27071153"
    citation: "Kuhle J, Barro C, Andreasson U, et al. Comparison of three analytical platforms for quantification of the neurofilament light chain in blood samples: ELISA, electrochemiluminescence immunoassay and Simoa. Clin Chem Lab Med. 2016;54(10):1655-1661."
    claim: "Analytical sensitivity for serum/CSF NfL quantification was 0.62 pg/mL for the Simoa (single-molecule array) digital immunoassay versus 78.0 pg/mL for conventional ELISA and 15.6 pg/mL for an electrochemiluminescence assay; only the more sensitive platforms reliably distinguished disease-related differences in blood NfL."
    verified: true
    verified_by: research-reviewer (claude@ogdb.com)

  - id: "doi:10.1001/jamaneurol.2021.3671"
    type: doi
    doi: "10.1001/jamaneurol.2021.3671"
    pmid: "34661615"
    citation: "Benedet AL, Milà-Alomà M, Vrillon A, et al. Differences Between Plasma and Cerebrospinal Fluid Glial Fibrillary Acidic Protein Levels Across the Alzheimer Disease Continuum. JAMA Neurol. 2021;78(12):1471-1483."
    claim: "Across three independent cohorts (TRIAD, ALFA+, BioCogBank Paris Lariboisière) spanning preclinical to symptomatic AD, plasma GFAP (measured with a Simoa assay) tracked the AD continuum and was associated with amyloid-PET/CSF Abeta42/40 status; GFAP is described as a marker of reactive astrogliosis elevated in both CSF and blood of individuals with AD, with plasma and CSF GFAP showing distinct performance characteristics across the continuum."
    verified: true
    verified_by: research-reviewer (claude@ogdb.com)

  - id: "doi:10.1038/s41398-020-01137-1"
    type: doi
    doi: "10.1038/s41398-020-01137-1"
    pmid: "33431793"
    citation: "Chatterjee P, Pedrini S, Stoops E, et al. Plasma glial fibrillary acidic protein is elevated in cognitively normal older adults at risk of Alzheimer's disease. Transl Psychiatry. 2021;11(1):27."
    claim: "In cognitively normal older adults (65-90y) classified by amyloid-PET, plasma GFAP (Simoa) was significantly higher in amyloid-positive (n=33) versus amyloid-negative (n=63) participants (P<0.00001); adding plasma GFAP to a covariate-only model raised AUC for discriminating amyloid status from 0.78 to 0.91, indicating astrocytic activation begins at the presymptomatic stage of AD."
    verified: true
    verified_by: research-reviewer (claude@ogdb.com)

  - id: "doi:10.1093/brain/awab223"
    type: doi
    doi: "10.1093/brain/awab223"
    pmid: "34259835"
    citation: "Pereira JB, Janelidze S, Smith R, et al. Plasma GFAP is an early marker of amyloid-beta but not tau pathology in Alzheimer's disease. Brain. 2021;144(11):3505-3516."
    claim: "In the Swedish BioFINDER-2 study (n=504), plasma GFAP concentration was significantly increased in all amyloid-beta-PET-positive groups (P<0.01) and predicted amyloid-PET positivity with AUC 0.76, outperforming CSF GFAP (AUC 0.69) and other glial markers including CSF YKL-40 (AUC 0.64) and soluble TREM2 (AUC 0.71); plasma GFAP tracked amyloid pathology but not tau pathology independent of amyloid."
    verified: true
    verified_by: research-reviewer (claude@ogdb.com)

  - id: "doi:10.1016/j.jalz.2015.05.012"
    type: doi
    doi: "10.1016/j.jalz.2015.05.012"
    pmid: "26092348"
    citation: "De Vos A, Jacobs D, Struyfs H, et al. C-terminal neurogranin is increased in cerebrospinal fluid but unchanged in plasma in Alzheimer's disease. Alzheimers Dement. 2015;11(12):1461-1469."
    claim: "In paired CSF/plasma samples from controls (n=29) and MCI/AD-dementia patients (n=59) measured with a C-terminal-neurogranin ELISA, CSF neurogranin was significantly increased in AD and correlated with CSF tau, but plasma neurogranin showed no difference between controls and AD, and CSF and plasma neurogranin did not correlate with one another."
    verified: true
    verified_by: research-reviewer (claude@ogdb.com)

  - id: "doi:10.1016/j.biopsych.2010.08.025"
    type: doi
    doi: "10.1016/j.biopsych.2010.08.025"
    pmid: "21035623"
    citation: "Craig-Schapiro R, Perrin RJ, Roe CM, et al. YKL-40: a novel prognostic fluid biomarker for preclinical Alzheimer's disease. Biol Psychiatry. 2010;68(10):903-912."
    claim: "CSF YKL-40 (identified by CSF proteomics, confirmed by ELISA in a discovery cohort n=47 and validation cohort n=292) was higher in very-mild/mild AD dementia (CDR 0.5/1) than controls (CDR 0), and the CSF YKL-40/Abeta42 ratio predicted conversion to cognitive impairment as well as CSF tau/Abeta42; in paired plasma samples (n=237) mean plasma YKL-40 was also higher in CDR 0.5/1 versus CDR 0 and correlated with CSF levels, though the primary discovery and validation was CSF-based."
    verified: true
    verified_by: research-reviewer (claude@ogdb.com)

  - id: "doi:10.1016/j.neurobiolaging.2020.07.009"
    type: doi
    doi: "10.1016/j.neurobiolaging.2020.07.009"
    pmid: "32920471"
    citation: "Vergallo A, Lista S, Lemercier P, et al; INSIGHT-preAD study group and Alzheimer Precision Medicine Initiative. Association of plasma YKL-40 with brain amyloid-beta levels, memory performance, and sex in subjective memory complainers. Neurobiol Aging. 2020;96:22-32."
    claim: "In a single monocentric longitudinal cohort of cognitively healthy individuals at risk of AD (INSIGHT-preAD), plasma YKL-40 increased with age, was higher in men than women, and was positively associated with memory performance but negatively associated with brain amyloid-beta deposition -- an inverse relationship to what would be expected of a marker of pathological progression, illustrating that plasma YKL-40 findings in AD-risk cohorts are inconsistent/study-specific rather than a validated directional biomarker."
    verified: true
    verified_by: research-reviewer (claude@ogdb.com)

  - id: "doi:10.1016/S1474-4422(16)00070-3"
    type: doi
    doi: "10.1016/S1474-4422(16)00070-3"
    pmid: "27068280"
    citation: "Olsson B, Lautner R, Andreasson U, et al. CSF and blood biomarkers for the diagnosis of Alzheimer's disease: a systematic review and meta-analysis. Lancet Neurol. 2016;15(7):673-684."
    claim: "In a systematic review/meta-analysis of 231 articles (15,699 AD patients, 13,018 controls) covering 15 CSF and blood biomarkers (including glial activation markers YKL-40 and GFAP), CSF NFL (ratio 2.35) and plasma T-tau (ratio 1.95) had large effect sizes distinguishing AD from controls, whereas CSF YKL-40 (along with NSE, VLP-1, HFABP) had only a moderate effect size (average ratio 1.28-1.47); the core recommendation was that CSF T-tau, P-tau, Abeta42 and NFL -- not blood glial markers -- were consistent enough for clinical/research use as of this 2016 analysis."
    verified: true
    verified_by: research-reviewer (claude@ogdb.com)

  - id: "doi:10.1016/j.jalz.2018.08.006"
    type: doi
    doi: "10.1016/j.jalz.2018.08.006"
    pmid: "30321501"
    citation: "Sandelius A, Portelius E, Kallen A, et al. Elevated CSF GAP-43 is Alzheimer's disease specific and associated with tau and amyloid pathology. Alzheimers Dement. 2019;15(1):55-64."
    claim: "CSF GAP-43 (growth-associated protein 43 / neuromodulin) was elevated in AD dementia and MCI-due-to-AD compared to controls and to non-AD neurodegenerative disease, was AD-specific among the dementias tested, and correlated with CSF tau and amyloid pathology; this study measured GAP-43 in CSF only, not blood/plasma."
    verified: true
    verified_by: research-reviewer (claude@ogdb.com)

  - id: "doi:10.1186/1750-1326-9-53"
    type: doi
    doi: "10.1186/1750-1326-9-53"
    pmid: "25418885"
    citation: "Brinkmalm A, Brinkmalm G, Honer WG, et al. SNAP-25 is a promising novel cerebrospinal fluid biomarker for synapse degeneration in Alzheimer's disease. Mol Neurodegener. 2014;9:53."
    claim: "SNAP-25 was identified and quantified as a synaptic-degeneration biomarker in CSF and proposed as a promising marker of synapse loss in AD; this original characterization was performed in CSF, not blood/plasma."
    verified: true
    verified_by: research-reviewer (claude@ogdb.com)

  - id: "doi:10.1038/nbt.1641"
    type: doi
    doi: "10.1038/nbt.1641"
    pmid: "20495550"
    citation: "Rissin DM, Kan CW, Campbell TG, et al. Single-molecule enzyme-linked immunosorbent assay detects serum proteins at subfemtomolar concentrations. Nat Biotechnol. 2010;28(6):595-599."
    claim: "The Simoa (single-molecule array) digital ELISA technology can quantify serum/plasma proteins down to subfemtomolar concentrations, roughly 1,000-fold more sensitive than conventional ELISA -- the technology class subsequently used to make plasma NfL and plasma GFAP measurable in AD cohorts (see Kuhle et al. 2016; Mattsson et al. 2017; Benedet et al. 2021; Chatterjee et al. 2021; Pereira et al. 2021)."
    verified: true
    verified_by: research-reviewer (claude@ogdb.com)

  - id: "doi:10.1074/mcp.r200007-mcp200"
    type: doi
    doi: "10.1074/mcp.r200007-mcp200"
    pmid: "12488461"
    citation: "Anderson NL, Anderson NG. The human plasma proteome: history, character, and diagnostic prospects. Mol Cell Proteomics. 2002;1(11):845-867."
    claim: "The human plasma proteome spans an estimated ~10 orders of magnitude in protein concentration, dominated by a small number of very-high-abundance proteins (e.g., albumin, immunoglobulins); this dynamic range is the central analytical obstacle to detecting low-abundance, tissue-leakage signaling proteins (the class NfL and GFAP belong to) by direct, untargeted plasma proteomic methods without depletion or enrichment."
    verified: true
    verified_by: research-reviewer (claude@ogdb.com)

  - id: "url:https://www.uniprot.org/uniprotkb/P07196/entry"
    type: url
    citation: "UniProtKB entry P07196 (NFL_HUMAN), Neurofilament light polypeptide, Homo sapiens, gene NEFL. Accessed 2026-09-15."
    claim: "Confirms canonical human UniProt accession P07196, entry name NFL_HUMAN, gene name NEFL, for neurofilament light polypeptide."
    verified: true
    verified_by: research-reviewer (claude@ogdb.com)

  - id: "url:https://www.uniprot.org/uniprotkb/P14136/entry"
    type: url
    citation: "UniProtKB entry P14136 (GFAP_HUMAN), Glial fibrillary acidic protein, Homo sapiens, gene GFAP. Accessed 2026-09-15."
    claim: "Confirms canonical human UniProt accession P14136, entry name GFAP_HUMAN, gene name GFAP, for glial fibrillary acidic protein."
    verified: true
    verified_by: research-reviewer (claude@ogdb.com)

  - id: "url:https://www.uniprot.org/uniprotkb/Q92686/entry"
    type: url
    citation: "UniProtKB entry Q92686 (NEUG_HUMAN), Neurogranin, Homo sapiens, gene NRGN. Accessed 2026-09-15."
    claim: "Confirms canonical human UniProt accession Q92686 for neurogranin, gene name NRGN. NOTE: the UniProt entry name is NEUG_HUMAN, not NRGN_HUMAN."
    verified: true
    verified_by: research-reviewer (claude@ogdb.com)

  - id: "url:https://www.uniprot.org/uniprotkb/P36222/entry"
    type: url
    citation: "UniProtKB entry P36222 (CH3L1_HUMAN), Chitinase-3-like protein 1 (YKL-40), Homo sapiens, gene CHI3L1. Accessed 2026-09-15."
    claim: "Confirms canonical human UniProt accession P36222 for chitinase-3-like protein 1 (YKL-40), gene name CHI3L1. NOTE: the UniProt entry name is CH3L1_HUMAN, not CHI3L1_HUMAN."
    verified: true
    verified_by: research-reviewer (claude@ogdb.com)

  - id: "url:https://www.uniprot.org/uniprotkb/P17677/entry"
    type: url
    citation: "UniProtKB entry P17677 (NEUM_HUMAN), Neuromodulin (GAP-43), Homo sapiens, gene GAP43. Accessed 2026-09-15."
    claim: "Confirms canonical human UniProt accession P17677 for neuromodulin/GAP-43, gene name GAP43. NOTE: the UniProt entry name is NEUM_HUMAN, not GAP43_HUMAN."
    verified: true
    verified_by: research-reviewer (claude@ogdb.com)

  - id: "url:https://www.uniprot.org/uniprotkb/P60880/entry"
    type: url
    citation: "UniProtKB entry P60880 (SNP25_HUMAN), Synaptosomal-associated protein 25 (SNAP-25), Homo sapiens, gene SNAP25. Accessed 2026-09-15."
    claim: "Confirms canonical human UniProt accession P60880 for SNAP-25, gene name SNAP25. NOTE: the UniProt entry name is SNP25_HUMAN, not SNAP25_HUMAN."
    verified: true
    verified_by: research-reviewer (claude@ogdb.com)
---

# Research: Glial, axonal, and synaptic neurodegeneration protein biomarkers of Alzheimer's disease in blood plasma

## Summary

Of the four core proteins in scope, two have strong, repeatedly replicated **plasma** evidence in AD and are appropriate to treat as genuine plasma neurodegeneration/glial-activation markers: **neurofilament light chain (NEFL/NfL)**, the canonical non-specific axonal-injury marker, elevated in plasma across the AD continuum and predictive of neurodegeneration and clinical progression ([doi:10.1001/jamaneurol.2016.6117]; [doi:10.1038/s41591-018-0304-3]); and **GFAP**, the canonical astrocytic-reactivity marker, elevated in plasma from the preclinical (amyloid-positive, cognitively normal) stage onward and, in at least one large cohort, a *better* predictor of amyloid-PET positivity than CSF GFAP or CSF YKL-40 ([doi:10.1001/jamaneurol.2021.3671]; [doi:10.1038/s41398-020-01137-1]; [doi:10.1093/brain/awab223]). The other two are markedly weaker in blood: **neurogranin (NRGN)**, a well-established CSF synaptic-loss marker, was explicitly reported as *unchanged* in plasma in a direct paired CSF/plasma comparison ([doi:10.1016/j.jalz.2015.05.012]) and there is no subsequent large positive plasma replication located in this search. **CHI3L1/YKL-40**, a glial/neuroinflammation marker, is predominantly a CSF finding ([doi:10.1016/j.biopsych.2010.08.025]; [doi:10.1016/S1474-4422(16)00070-3]); plasma YKL-40 data exist but are inconsistent in direction and effect across studies, including one report of an *inverse* (protective-looking) association with brain amyloid ([doi:10.1016/j.neurobiolaging.2020.07.009]) and only moderate performance versus GFAP as a plasma glial marker in a head-to-head comparison ([doi:10.1093/brain/awab223]).

For the QC dynamic-range/MS-detectability question: plasma NEFL and GFAP concentrations in these AD studies sit in the **ng/L (i.e., low pg/mL) range** (e.g., plasma NfL means of 34.7-51.0 ng/L in ADNI; [doi:10.1001/jamaneurol.2016.6117]), and reaching that sensitivity in blood required moving to digital/single-molecule immunoassay (Simoa) technology, which is roughly 1,000-fold more sensitive than conventional ELISA (Simoa analytical sensitivity 0.62 pg/mL vs 78.0 pg/mL for ELISA; [doi:10.1515/cclm-2015-1195]; [doi:10.1038/nbt.1641]). Combined with the ~10-order-of-magnitude dynamic range of the plasma proteome dominated by a handful of very-high-abundance proteins ([doi:10.1074/mcp.r200007-mcp200]), this means GFAP and NfL sit deep in the low-abundance tail of the plasma proteome, and untargeted, non-enriched plasma DIA-MS runs at typical proteomic depth commonly do **not** reach the sensitivity of Simoa; whether they are detected in an EV-enriched plasma DIA-MS workflow specifically was not established by any source found in this search and should be treated as an open, empirical question for the QC plot rather than an assumed positive.

GAP43 and SNAP25 are additional synaptic-loss markers with solid CSF evidence ([doi:10.1016/j.jalz.2018.08.006]; [doi:10.1186/1750-1326-9-53]) but no blood/plasma evidence was located in this search; they should be labeled CSF-only if referenced.

## Detailed findings

### 1. NEFL / neurofilament light chain — STRONG plasma evidence (non-specific axonal-injury marker)

- **Identity (verified against UniProt):** gene `NEFL`; human UniProt accession **P07196**; entry name **NFL_HUMAN** (`sp|P07196|NFL_HUMAN`); recommended name "Neurofilament light polypeptide" ([url:https://www.uniprot.org/uniprotkb/P07196/entry]).
- **Evidence context: PLASMA/BLOOD (strong, multi-cohort).** In the ADNI cohort (n=570: 193 controls, 197 MCI, 180 AD dementia), plasma NfL measured with an ultrasensitive assay was significantly elevated in MCI (mean 42.8 ng/L) and AD dementia (mean 51.0 ng/L) versus controls (mean 34.7 ng/L, P<.001), gave AUC 0.87 for AD dementia vs controls (comparable to established CSF biomarkers), correlated with CSF NfL (rho=0.59), and tracked cognitive decline, AD-pattern atrophy, and hypometabolism longitudinally ([doi:10.1001/jamaneurol.2016.6117]).
- In the Dominantly Inherited Alzheimer Network (autosomal-dominant/familial AD), serum NfL (n=405 cross-sectional, n=196 longitudinal) was elevated at presymptomatic stages and correlated with CSF NfL; critically, the *rate of change* of serum NfL discriminated mutation carriers from non-carriers ~16.2 years before estimated symptom onset, versus only ~6.8 years using a single cross-sectional NfL level — i.e., serial plasma/serum NfL adds substantial lead time over a single measurement ([doi:10.1038/s41591-018-0304-3]).
- NfL is explicitly **non-specific**: it is a marker of axonal injury/neurodegeneration seen across many neurological diseases (both source studies frame it this way — DIAN is a familial-AD-specific cohort but NfL itself is described as "a promising fluid biomarker of disease progression for various cerebral proteopathies," and ADNI shows it also rises in MCI, not just dementia) ([doi:10.1038/s41591-018-0304-3]; [doi:10.1001/jamaneurol.2016.6117]). It should be interpreted as a general neurodegeneration ("N") marker, not an AD-specific pathology marker.
- **MS-detectability:** Reaching disease-relevant blood NfL concentrations required moving from conventional ELISA (analytical sensitivity 78.0 pg/mL) to Simoa digital ELISA (0.62 pg/mL); an electrochemiluminescence platform was intermediate (15.6 pg/mL) ([doi:10.1515/cclm-2015-1195]). Given plasma NfL concentrations in AD cohorts of ~35-51 ng/L (= pg/mL range), NfL is an ultralow-abundance plasma protein by mass-spectrometry standards.

### 2. GFAP / glial fibrillary acidic protein — STRONG plasma evidence (astrocytic-reactivity marker)

- **Identity (verified against UniProt):** gene `GFAP`; human UniProt accession **P14136**; entry name **GFAP_HUMAN** (`sp|P14136|GFAP_HUMAN`); recommended name "Glial fibrillary acidic protein" ([url:https://www.uniprot.org/uniprotkb/P14136/entry]).
- **Evidence context: PLASMA/BLOOD (strong, multi-cohort, includes head-to-head vs CSF).** Across three independent cohorts (TRIAD [Montreal], ALFA+ [Barcelona], BioCogBank Paris Lariboisière), plasma GFAP (Simoa assay) tracked the entire AD continuum from preclinical to dementia and was associated with amyloid pathology; the study explicitly frames GFAP as "a marker of reactive astrogliosis that increases in the cerebrospinal fluid (CSF) and blood of individuals with Alzheimer disease" and evaluates whether plasma performance matches CSF performance across stages ([doi:10.1001/jamaneurol.2021.3671]).
- In cognitively **normal** older adults classified only by amyloid-PET status, plasma GFAP was significantly higher in amyloid-positive (n=33) than amyloid-negative (n=63) participants (P<0.00001); adding plasma GFAP to a demographic/APOE covariate model raised discrimination of amyloid status from AUC 0.78 to 0.91 — i.e., astrocytic activation detectable in plasma precedes cognitive symptoms ([doi:10.1038/s41398-020-01137-1]).
- In the Swedish BioFINDER-2 study (n=504 across the amyloid/tau/cognition spectrum), plasma GFAP was significantly elevated in every amyloid-PET-positive group (P<0.01) and predicted amyloid-PET positivity with **AUC 0.76**, which was *higher* than CSF GFAP (AUC 0.69) and higher than CSF YKL-40 (AUC 0.64) or soluble TREM2 (AUC 0.71) measured in the same cohort. Plasma GFAP tracked amyloid pathology specifically, with associations to tau-PET becoming non-significant after controlling for amyloid-PET ([doi:10.1093/brain/awab223]). This is a directly relevant finding for prioritizing GFAP over YKL-40 as the plasma glial marker of choice.
- **MS-detectability:** All three plasma GFAP studies cited here used the Simoa digital immunoassay platform, not conventional immunoassay or MS, implying GFAP — like NfL — sits at very low (sub-ELISA-sensitivity) concentrations in plasma ([doi:10.1001/jamaneurol.2021.3671]; [doi:10.1038/s41398-020-01137-1]; [doi:10.1093/brain/awab223]; general Simoa sensitivity context: [doi:10.1038/nbt.1641]).

### 3. NRGN / neurogranin — WEAK/inconsistent plasma evidence; strong CSF-only marker

- **Identity (verified against UniProt):** gene `NRGN`; human UniProt accession **Q92686**; entry name **NEUG_HUMAN** (`sp|Q92686|NEUG_HUMAN`) — note this differs from the naive expectation of "NRGN_HUMAN"; recommended name "Neurogranin" ([url:https://www.uniprot.org/uniprotkb/Q92686/entry]).
- **Evidence context: CSF strong / PLASMA weak-to-absent.** In a direct, paired CSF-and-plasma comparison (controls n=29; MCI/AD-dementia n=59, using a C-terminal-neurogranin ELISA), **CSF neurogranin was significantly increased in AD** and correlated with CSF tau, while **plasma neurogranin showed no difference between controls and AD**, and CSF and plasma neurogranin levels did not correlate with each other, arguing against plasma neurogranin simply reflecting the CSF signal at lower amplitude ([doi:10.1016/j.jalz.2015.05.012]). The paper's own framing states neurogranin "may be a promising cerebrospinal fluid (CSF) biomarker but has not yet been evaluated as a plasma biomarker" prior to their study, and their conclusion "strengthens the potential of neurogranin as an AD **CSF** biomarker" specifically (not a blood biomarker).
- No subsequent large positive plasma-neurogranin replication was located in this targeted search; treat plasma/blood NRGN as **not established** and CSF-dominant. If NRGN is included in a plasma DIA-MS panel, the QC narrative should state explicitly that its inclusion tests a marker with a documented negative plasma finding, not an established positive one.

### 4. CHI3L1 / YKL-40 — CSF-dominant evidence; plasma data exist but are inconsistent

- **Identity (verified against UniProt):** gene `CHI3L1`; human UniProt accession **P36222**; entry name **CH3L1_HUMAN** (`sp|P36222|CH3L1_HUMAN`) — note this differs from the naive expectation of "CHI3L1_HUMAN"; recommended name "Chitinase-3-like protein 1"; alternative name "YKL-40" ([url:https://www.uniprot.org/uniprotkb/P36222/entry]).
- **Evidence context: CSF strong / PLASMA present but weaker and inconsistent.** The original characterization study used CSF proteomics for discovery and CSF ELISA for validation (discovery n=47, validation n=292); CSF YKL-40 was higher in very-mild/mild AD dementia (CDR 0.5/1) than controls, and CSF YKL-40/Abeta42 predicted conversion to cognitive impairment about as well as CSF tau/Abeta42 or p-tau181/Abeta42. The same study *did* also measure paired plasma samples (n=237) and reported mean plasma YKL-40 was higher in CDR 0.5/1 versus CDR 0 and correlated with CSF levels — so there is direct-comparison plasma support, but it is secondary to (and smaller/less emphasized than) the CSF finding, from a single study using conventional (non-ultrasensitive) ELISA ([doi:10.1016/j.biopsych.2010.08.025]).
- A 2016 systematic review/meta-analysis of 231 studies (15,699 AD patients, 13,018 controls) covering 15 CSF/blood biomarkers including the glial-activation markers YKL-40 and GFAP found **CSF YKL-40 had only a moderate effect size** (average ratio 1.28-1.47, in the same tier as NSE/VLP-1/HFABP) distinguishing AD from controls — clearly weaker than CSF T-tau/P-tau/Abeta42 or CSF NFL (ratio 2.35) — and did not identify blood/plasma YKL-40 among the markers with a consistent, strong AD-associated signal; the review's overall recommendation was to rely on the CSF core biomarkers plus CSF NFL, not blood glial markers ([doi:10.1016/S1474-4422(16)00070-3]).
- A later single-cohort study (INSIGHT-preAD, cognitively healthy individuals at risk of AD) found plasma YKL-40 increased with age and was higher in men, and was *negatively* associated with brain amyloid-beta deposition and *positively* associated with memory performance — i.e., in this cohort higher plasma YKL-40 tracked with *better*, not worse, AD-pathology status, the opposite direction from what a disease-progression marker would be expected to show. The authors interpret this as a possible early compensatory/protective glial response rather than treating it as a straightforward positive AD biomarker ([doi:10.1016/j.neurobiolaging.2020.07.009]).
- In the BioFINDER-2 head-to-head comparison of glial markers for amyloid-PET prediction, **CSF YKL-40 gave the weakest AUC (0.64)** of the glial/neuroinflammation markers tested (vs plasma GFAP AUC 0.76, CSF GFAP AUC 0.69, soluble TREM2 AUC 0.71); plasma YKL-40 was not separately reported as a top-performing marker in that comparison ([doi:10.1093/brain/awab223]).
- **Net assessment:** CHI3L1/YKL-40 should be labeled **CSF-dominant, with plasma evidence that is real but inconsistent in magnitude and even direction across studies** — materially weaker and less reproducible than plasma GFAP.

### 5. Optional CSF-only synaptic markers: GAP43 and SNAP25 (no blood/plasma evidence located)

- **GAP43 / neuromodulin.** Gene `GAP43`; human UniProt accession **P17677**; entry name **NEUM_HUMAN** (`sp|P17677|NEUM_HUMAN`) ([url:https://www.uniprot.org/uniprotkb/P17677/entry]). CSF GAP-43 was elevated in AD dementia and MCI-due-to-AD versus controls and versus non-AD neurodegenerative disease (i.e., AD-specific among the dementias tested), and correlated with CSF tau and amyloid pathology. This characterization was performed **in CSF only**; no blood/plasma measurement was part of this study ([doi:10.1016/j.jalz.2018.08.006]).
- **SNAP25.** Gene `SNAP25`; human UniProt accession **P60880**; entry name **SNP25_HUMAN** (`sp|P60880|SNP25_HUMAN`) ([url:https://www.uniprot.org/uniprotkb/P60880/entry]). SNAP-25 was identified and validated as a synapse-degeneration biomarker **in CSF**, described as "a promising novel cerebrospinal fluid biomarker for synapse degeneration in Alzheimer's disease." No blood/plasma measurement was part of the original characterization ([doi:10.1186/1750-1326-9-53]).
- Neither protein has a plasma AD biomarker study identified in this targeted search; if mentioned in the manuscript at all, they should be labeled explicitly as **CSF-only / not established in blood**.

### 6. MS-detectability context for the QC dynamic-range plot

- The human plasma proteome spans an estimated **~10 orders of magnitude** in protein concentration, dominated by a small number of very-high-abundance proteins (e.g., albumin, immunoglobulins); this dynamic range is the classic analytical barrier to detecting low-abundance, tissue-leakage signaling proteins in unfractionated plasma by direct proteomic methods ([doi:10.1074/mcp.r200007-mcp200]).
- NfL and GFAP are quantified in the AD literature almost exclusively via **Simoa (single-molecule array) digital immunoassay**, a technology built specifically to reach subfemtomolar/low-picomolar sensitivity — roughly 1,000-fold more sensitive than conventional ELISA ([doi:10.1038/nbt.1641]; head-to-head sensitivity numbers: Simoa 0.62 pg/mL vs ELISA 78.0 pg/mL vs electrochemiluminescence 15.6 pg/mL, [doi:10.1515/cclm-2015-1195]). Reported plasma NfL concentrations in AD cohorts are in the **ng/L (low pg/mL) range** (34.7-51.0 ng/L in ADNI; [doi:10.1001/jamaneurol.2016.6117]).
- **Implication for this DIA-MS study:** GFAP and NfL are textbook examples of ultralow-abundance plasma proteins that were historically Simoa/immunoassay territory precisely because standard proteomic workflows lack the dynamic range to reach them in neat, unfractionated plasma. None of the sources reviewed here report untargeted plasma (or plasma-EV) DIA-MS detection of GFAP or NfL, so **whether this study's EV-enriched plasma DIA-MS workflow detects them at all is an open empirical question** that the dynamic-range QC plot itself should answer — it is not something the literature can currently confirm or rule out either way. NRGN and CHI3L1/YKL-40 are comparatively higher-abundance/less specialized-assay-dependent in the sources reviewed (CHI3L1/YKL-40 has been measured by conventional ELISA in plasma; [doi:10.1016/j.biopsych.2010.08.025]), but this reflects assay history, not confirmed MS detectability.

## Methods / sources consulted

- **UniProt (REST API, `rest.uniprot.org/uniprotkb/<accession>.txt`)**: fetched and read the full flat-file record for each of P07196, P14136, Q92686, P36222 directly (confirming entry name, accession, gene name, organism = Homo sapiens) on 2026-09-15; searched by gene symbol (`gene:GAP43`, `gene:SNAP25`, human, reviewed/Swiss-Prot only) via `rest.uniprot.org/uniprotkb/search` to resolve P17677 and P60880. All six accessions and entry names in this document were read directly from UniProt, not recalled from memory — two (NRGN, CHI3L1) turned out to differ from the naively expected entry-name pattern (`NEUG_HUMAN`, `CH3L1_HUMAN` rather than `NRGN_HUMAN`, `CHI3L1_HUMAN`), which is exactly the kind of detail this verification step is meant to catch.
- **PubMed / NCBI E-utilities** (`esearch`, `esummary`, `efetch` against `eutils.ncbi.nlm.nih.gov`): used to locate candidate papers by targeted author/topic queries, then confirmed each with `esummary` (title, journal, year, DOI) and `efetch` (full abstract) before citing it, reading enough of the abstract in every case to confirm the specific numeric/directional claim attributed to it (not just title/topic match). No citation in this document was accepted on title match alone.
- **Search strategy**: targeted PubMed queries per protein (e.g., "plasma neurofilament light neurodegeneration Alzheimer", "plasma glial fibrillary acidic protein cognitively normal", "C-terminal neurogranin cerebrospinal fluid plasma Alzheimer", "plasma YKL-40 CHI3L1 Alzheimer") rather than a broad literature survey, per the bounded scope of this finding. A general CSF-vs-blood meta-analysis (Olsson et al. 2016, Lancet Neurology) was added for aggregate/comparative context across all four core markers.
- **Not done / limitations**: did not access full-text PDFs (abstracts + structured PubMed metadata only); did not perform a systematic/exhaustive literature search (e.g., no formal search-string documentation, no citation-count-based recency check for CHI3L1/YKL-40 beyond the papers found); did not locate a study directly reporting untargeted plasma or plasma-EV DIA-MS peptide-level detection of GFAP or NfL — this remains an explicit gap, flagged above rather than papered over.

## References
Rendered from `references` (see YAML frontmatter for structured entries; all currently `verified: false`, pending research-reviewer sign-off):

1. Mattsson N, Andreasson U, Zetterberg H, Blennow K; Alzheimer's Disease Neuroimaging Initiative. Association of Plasma Neurofilament Light With Neurodegeneration in Patients With Alzheimer Disease. *JAMA Neurol.* 2017;74(5):557-566. doi:10.1001/jamaneurol.2016.6117; PMID:28346578.
2. Preische O, Schultz SA, Apel A, et al; Dominantly Inherited Alzheimer Network. Serum neurofilament dynamics predicts neurodegeneration and clinical progression in presymptomatic Alzheimer's disease. *Nat Med.* 2019;25(2):277-283. doi:10.1038/s41591-018-0304-3; PMID:30664784.
3. Kuhle J, Barro C, Andreasson U, et al. Comparison of three analytical platforms for quantification of the neurofilament light chain in blood samples: ELISA, electrochemiluminescence immunoassay and Simoa. *Clin Chem Lab Med.* 2016;54(10):1655-1661. doi:10.1515/cclm-2015-1195; PMID:27071153.
4. Benedet AL, Milà-Alomà M, Vrillon A, et al. Differences Between Plasma and Cerebrospinal Fluid Glial Fibrillary Acidic Protein Levels Across the Alzheimer Disease Continuum. *JAMA Neurol.* 2021;78(12):1471-1483. doi:10.1001/jamaneurol.2021.3671; PMID:34661615.
5. Chatterjee P, Pedrini S, Stoops E, et al. Plasma glial fibrillary acidic protein is elevated in cognitively normal older adults at risk of Alzheimer's disease. *Transl Psychiatry.* 2021;11(1):27. doi:10.1038/s41398-020-01137-1; PMID:33431793.
6. Pereira JB, Janelidze S, Smith R, et al. Plasma GFAP is an early marker of amyloid-beta but not tau pathology in Alzheimer's disease. *Brain.* 2021;144(11):3505-3516. doi:10.1093/brain/awab223; PMID:34259835.
7. De Vos A, Jacobs D, Struyfs H, et al. C-terminal neurogranin is increased in cerebrospinal fluid but unchanged in plasma in Alzheimer's disease. *Alzheimers Dement.* 2015;11(12):1461-1469. doi:10.1016/j.jalz.2015.05.012; PMID:26092348.
8. Craig-Schapiro R, Perrin RJ, Roe CM, et al. YKL-40: a novel prognostic fluid biomarker for preclinical Alzheimer's disease. *Biol Psychiatry.* 2010;68(10):903-912. doi:10.1016/j.biopsych.2010.08.025; PMID:21035623.
9. Vergallo A, Lista S, Lemercier P, et al; INSIGHT-preAD study group and Alzheimer Precision Medicine Initiative. Association of plasma YKL-40 with brain amyloid-beta levels, memory performance, and sex in subjective memory complainers. *Neurobiol Aging.* 2020;96:22-32. doi:10.1016/j.neurobiolaging.2020.07.009; PMID:32920471.
10. Olsson B, Lautner R, Andreasson U, et al. CSF and blood biomarkers for the diagnosis of Alzheimer's disease: a systematic review and meta-analysis. *Lancet Neurol.* 2016;15(7):673-684. doi:10.1016/S1474-4422(16)00070-3; PMID:27068280.
11. Sandelius A, Portelius E, Kallen A, et al. Elevated CSF GAP-43 is Alzheimer's disease specific and associated with tau and amyloid pathology. *Alzheimers Dement.* 2019;15(1):55-64. doi:10.1016/j.jalz.2018.08.006; PMID:30321501.
12. Brinkmalm A, Brinkmalm G, Honer WG, et al. SNAP-25 is a promising novel cerebrospinal fluid biomarker for synapse degeneration in Alzheimer's disease. *Mol Neurodegener.* 2014;9:53. doi:10.1186/1750-1326-9-53; PMID:25418885.
13. Rissin DM, Kan CW, Campbell TG, et al. Single-molecule enzyme-linked immunosorbent assay detects serum proteins at subfemtomolar concentrations. *Nat Biotechnol.* 2010;28(6):595-599. doi:10.1038/nbt.1641; PMID:20495550.
14. Anderson NL, Anderson NG. The human plasma proteome: history, character, and diagnostic prospects. *Mol Cell Proteomics.* 2002;1(11):845-867. doi:10.1074/mcp.r200007-mcp200; PMID:12488461.
15. UniProtKB entry P07196 (NFL_HUMAN) — NEFL, Homo sapiens. https://www.uniprot.org/uniprotkb/P07196/entry (accessed 2026-09-15).
16. UniProtKB entry P14136 (GFAP_HUMAN) — GFAP, Homo sapiens. https://www.uniprot.org/uniprotkb/P14136/entry (accessed 2026-09-15).
17. UniProtKB entry Q92686 (NEUG_HUMAN) — NRGN, Homo sapiens. https://www.uniprot.org/uniprotkb/Q92686/entry (accessed 2026-09-15).
18. UniProtKB entry P36222 (CH3L1_HUMAN) — CHI3L1, Homo sapiens. https://www.uniprot.org/uniprotkb/P36222/entry (accessed 2026-09-15).
19. UniProtKB entry P17677 (NEUM_HUMAN) — GAP43, Homo sapiens. https://www.uniprot.org/uniprotkb/P17677/entry (accessed 2026-09-15).
20. UniProtKB entry P60880 (SNP25_HUMAN) — SNAP25, Homo sapiens. https://www.uniprot.org/uniprotkb/P60880/entry (accessed 2026-09-15).
