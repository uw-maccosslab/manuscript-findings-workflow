---
# Research-finding document — external knowledge, structured like a finding but for the literature/tools.
# Spec: doc 04.4. Stored in research/<slug>.md. A research finding without verified references is NOT accepted.

topic: "Core ATN amyloid (A) and tau (T) plasma protein biomarkers of Alzheimer's disease — APP-derived Abeta42/Abeta40 and MAPT-derived total-tau/p-tau181/p-tau217/p-tau231 — with UniProt identity verification and plasma DIA-MS detectability context for a plasma-EV study"
slug: ad-plasma-atn-core-biomarkers
type: research-publications

created: 2026-09-15
updated: 2026-09-15

status: reviewed
reviewed_by: research-reviewer (claude@ogdb.com)
reviewed_date: 2026-09-15

# Normalized domain-entity references this research is about (same scheme as findings; canonical IDs).
entities:
  - { db: uniprot, id: "P05067", label: "APP" }
  - { db: hgnc, id: "HGNC:620", label: "APP" }
  - { db: uniprot, id: "P10636", label: "MAPT" }
  - { db: hgnc, id: "HGNC:6893", label: "MAPT" }

# (required, non-empty) Every external claim traces to a reference here. Each is fact-checked by the
# research-reviewer: it must EXIST and SUPPORT the claim attributed to it. Hallucinated citations are
# the known failure mode this section + review exist to stop.
references:
  - id: "uniprot:P05067"
    type: url
    claim: "Human APP (gene symbol APP, HGNC:620) is UniProt accession P05067, entry name A4_HUMAN (canonical 770-aa sequence) — NOT 'APP_HUMAN'. Feature IDs for this protein in a bottom-up proteomics search database will read sp|P05067|A4_HUMAN. The Abeta42 biomarker peptide is chain/proteolytic-fragment residues 672-713 (UniProt feature id PRO_0000000092) and Abeta40 is residues 672-711 (PRO_0000000093) of this canonical sequence."
    verified: true
    verified_by: research-reviewer (claude@ogdb.com)

  - id: "uniprot:P10636"
    type: url
    claim: "Human MAPT (gene symbol MAPT, HGNC:6893) is UniProt accession P10636, entry name TAU_HUMAN (matches feature-id convention sp|P10636|TAU_HUMAN). The UniProt canonical/displayed sequence is isoform 'PNS-tau' (P10636-1, 758 aa) — not the 441-residue 2N4R 'Tau-F' isoform (P10636-8) that the AD biomarker literature conventionally uses for phospho-site numbering. Sequence alignment (performed directly against the UniProt FASTA records for P10636 and isoform P10636-8) shows Tau-441-numbering Thr181, Thr217, and Thr231 correspond to canonical UniProt residues Thr498, Thr534, and Thr548 respectively; all three are independently annotated as MOD_RES phosphothreonine sites in the live P10636 UniProt feature table, cross-validating the alignment-derived mapping."
    verified: true
    verified_by: research-reviewer (claude@ogdb.com)

  - id: "pmid:29420472"
    type: pmid
    claim: "Plasma amyloid-beta biomarkers derived from APP — the APP669-711/Abeta1-42 ratio and the Abeta1-40/Abeta1-42 ratio, plus their composite — predicted individual brain amyloid-PET status with high accuracy (AUC 96.7% discovery, 94.1% validation; ~90% accuracy vs PiB-PET). Measurement required immunoprecipitation coupled with mass spectrometry (IP-MS/MALDI-TOF), not untargeted/shotgun plasma proteomics."
    verified: true
    verified_by: research-reviewer (claude@ogdb.com)

  - id: "pmid:31371569"
    type: pmid
    claim: "A high-precision immunoprecipitation-and-LC-MS/MS assay for plasma Abeta42/Abeta40 (Washington University/C2N methodology) diagnosed brain amyloidosis with AUC 0.88 alone and AUC 0.94 combined with age and APOE4 status; a positive plasma Abeta42/40 despite negative baseline amyloid-PET conferred 15-fold greater risk of future PET conversion. The assay is explicitly an 'immunoprecipitation and liquid chromatography-mass spectrometry assay,' i.e. targeted immunoaffinity enrichment prior to MS, not untargeted shotgun/DIA."
    verified: true
    verified_by: research-reviewer (claude@ogdb.com)

  - id: "pmid:38595851"
    type: pmid
    claim: "Even with a purpose-built, high-precision LC-MS/MS assay for plasma Abeta42/40, the separation between amyloid-PET-positive and -negative groups was only ~18% (described as 'a relatively large difference... compared to most literature reports'), versus a typical ~8-11% separation reported for other published LC-MS/MS plasma Abeta42/40 assays; a 10% assay imprecision was shown elsewhere to produce a 26% misclassification rate for immunoassay-based Abeta42/40. This illustrates that the plasma Abeta42/40 signal is analytically a small fold-change, not a large-effect biomarker."
    verified: true
    verified_by: research-reviewer (claude@ogdb.com)

  - id: "pmid:16837572"
    type: pmid
    claim: "Beta-cleavage of APP occurs in early endosomes with routing of Abeta to multivesicular bodies, and a minor fraction of Abeta peptides is secreted in association with exosomes in HeLa/N2a cell culture; exosomal proteins were also found to accumulate in AD patient brain plaques. This is cell-culture and human brain-tissue evidence for an EV-association mechanism for APP/Abeta — it is NOT evidence that Abeta42/40 is a validated, quantifiable analyte in isolated human blood/plasma EVs; no plasma-EV-specific clinical biomarker performance is reported in this study."
    verified: true
    verified_by: research-reviewer (claude@ogdb.com)

  - id: "pmid:32333900"
    type: pmid
    claim: "An ultrasensitive blood immunoassay for plasma p-tau181 (MAPT phosphorylated at Thr181) distinguished Alzheimer's disease dementia from cognitively unimpaired older adults (AUC 90.2-98.2% across cohorts, n=1131 total) and from other neurodegenerative disorders (frontotemporal dementia AUC 82.8-100%, vascular dementia AUC 92.1%, PSP/CBS AUC 88.5%, PD/MSA AUC 81.9%), and predicted 1-year cognitive decline and hippocampal atrophy. Measured by ultrasensitive immunoassay (Simoa), in plasma, not mass spectrometry."
    verified: true
    verified_by: research-reviewer (claude@ogdb.com)

  - id: "pmid:32722745"
    type: pmid
    claim: "Plasma p-tau217 (MAPT phosphorylated at Thr217) discriminated neuropathologically- or clinically-defined Alzheimer disease from other neurodegenerative disorders across three independent cohorts (n=1402 total) with AUC 0.89-0.96, significantly higher accuracy than plasma p-tau181 and plasma neurofilament light chain (AUC 0.50-0.81, p<.001-.05), and performance not significantly different from CSF p-tau217/p-tau181 or tau-PET (AUC 0.90-0.99). In an autosomal-dominant AD kindred, plasma P-tau217 was elevated from ~25 years of age, ~20 years before estimated MCI onset in mutation carriers. Currently the single strongest plasma A/T-axis biomarker reported for AD discrimination."
    verified: true
    verified_by: research-reviewer (claude@ogdb.com)

  - id: "pmid:33585983"
    type: pmid
    claim: "Plasma p-tau231 (MAPT phosphorylated at Thr231), measured by an ultrasensitive Simoa immunoassay validated in 4 cohorts (n=588), identified AD and distinguished it from amyloid-negative cognitively unimpaired adults (AUC 0.92-0.94), from non-AD neurodegenerative disorders (AUC 0.93), and from amyloid-negative MCI (AUC 0.89). Plasma p-tau231 rose earlier than plasma p-tau181 and earlier than the standard amyloid-PET positivity threshold, and was elevated across the full Braak-stage spectrum including Braak 0-II, unlike plasma p-tau181 — positioning it as an early/incipient-pathology plasma marker."
    verified: true
    verified_by: research-reviewer (claude@ogdb.com)

  - id: "pmid:27694257"
    type: pmid
    claim: "Plasma total tau (unmodified MAPT) was associated with AD dementia status and with CSF tau/Abeta42 at the group level in a 1284-participant study (ADNI + BioFINDER), but the correlations were weak and differed between cohorts, with large overlap between AD and normal aging. The authors explicitly concluded these results 'do not support plasma tau as an AD biomarker in individual people' — i.e., of the tau-axis plasma markers, unmodified total tau is comparatively weak/nonspecific relative to the phospho-tau species."
    verified: true
    verified_by: research-reviewer (claude@ogdb.com)

  - id: "pmid:32725127"
    type: pmid
    claim: "Plasma p-tau217 and p-tau181 were quantified directly by mass spectrometry and were highly specific for amyloid-plaque pathology (AUROC 0.99 and 0.98 respectively in a 36-subject discovery cohort; 0.92 and 0.75 in a 92-subject validation cohort). Achieving this required that 'plasma tau was truncated' and the sample subsequently 'immunopurified against N-terminus and mid-domain tau antibodies' before MS — i.e., targeted immunoaffinity enrichment coupled to MS, not untargeted/global plasma DIA-MS."
    verified: true
    verified_by: research-reviewer (claude@ogdb.com)

  - id: "pmid:12488461"
    type: pmid
    claim: "The human plasma proteome has 'an extraordinary dynamic range in that more than 10 orders of magnitude in concentration separate albumin and the rarest proteins now measured clinically,' and the restricted dynamic range of conventional MS-based proteomic technology has historically limited plasma protein discovery to a small cataloged set — the dynamic-range problem that governs whether a low-abundance, non-plasma-resident (e.g., neuron-shed) protein/peptide is visible in an untargeted plasma digest without enrichment."
    verified: true
    verified_by: research-reviewer (claude@ogdb.com)

  - id: "pmid:27135364"
    type: pmid
    claim: "Single-run, non-depletion shotgun/LC-MS 'plasma proteome profiling' quantifies on the order of hundreds of plasma proteins per run (up to ~1,000 with peptide pre-fractionation), reproducibly covering the apolipoprotein family, inflammatory markers, and >40 FDA-approved biomarkers — i.e., analytes dominated by moderate-to-high-abundance plasma proteins, not brain-/neuron-derived low-abundance proteins."
    verified: true
    verified_by: research-reviewer (claude@ogdb.com)
---

# Research: Core ATN amyloid (A) and tau (T) plasma protein biomarkers of Alzheimer's disease (APP, MAPT)

## Summary

The plasma "A" and "T" axes of the ATN framework reduce, at the protein-identity level, to exactly two parent gene products: **APP** (UniProt `P05067`, entry name **`A4_HUMAN`** — not `APP_HUMAN`) supplies the amyloid-beta biomarker peptides **Abeta42** and the **Abeta42/Abeta40 ratio**; **MAPT** (UniProt `P10636`, entry name **`TAU_HUMAN`**) supplies **total tau** and the phosphorylated-tau species **p-tau181, p-tau217, and p-tau231** [uniprot:P05067; uniprot:P10636]. Across the published plasma literature, evidence strength is uneven: **plasma p-tau217 is currently the single strongest ATN plasma biomarker** (AUC up to 0.96-0.99 vs. other neurodegenerative disease and vs. CSF/PET) [pmid:32722745], **plasma p-tau181 is also strong and well-replicated** (AUC 0.90-0.98) [pmid:32333900], **plasma p-tau231 is a promising, earlier-rising marker** but less extensively replicated [pmid:33585983], and **plasma total tau is comparatively weak** — the primary study concluded it does not support individual-level AD diagnosis [pmid:27694257]. **Plasma Abeta42/Abeta40 is real and PET-concordant but analytically hard**: every study that reached strong performance used **immunoprecipitation coupled to mass spectrometry (IP-MS)**, not simple immunoassay or untargeted shotgun MS [pmid:29420472; pmid:31371569], and even the best LC-MS/MS assays report only an ~8-18% separation between amyloid-PET-positive and -negative groups [pmid:38595851]. Critically for this project's plasma-EV DIA-MS design: **none of the plasma p-tau or Abeta42/40 assays cited here — immunoassay or MS-based — used untargeted/global (non-enriched) plasma proteomics.** Every MS-based plasma quantification of these markers required immunoaffinity enrichment first [pmid:29420472; pmid:31371569; pmid:32725127], and the general plasma-proteome literature independently establishes why: the plasma proteome spans >10 orders of magnitude in abundance [pmid:12488461], and untargeted single-run LC-MS plasma profiling reproducibly reaches only a few hundred to ~1,000 proteins, dominated by moderate/high-abundance analytes [pmid:27135364] — several orders of magnitude short of where Abeta and phospho-tau peptides sit. No plasma-EV-specific (as opposed to whole/bulk plasma) clinical-performance literature for these exact ATN epitopes was located in this search; the closest mechanistic link found is cell-culture/brain-tissue evidence that Abeta can be exosome-associated [pmid:16837572], which is not evidence of a validated blood-EV biomarker.

## Detailed findings

### 1. Amyloid axis — parent protein APP

**Identity (verified against live UniProt record, 2026-09-15):**

| Field | Value | Source |
|---|---|---|
| Gene symbol | `APP` | UniProt GN line, ECO:0000312\|HGNC:HGNC:620 [uniprot:P05067] |
| UniProt accession (human) | `P05067` | [uniprot:P05067] |
| UniProt entry name | **`A4_HUMAN`** | [uniprot:P05067] |
| Canonical sequence length | 770 aa | [uniprot:P05067] |
| Recommended full name | Amyloid-beta precursor protein | [uniprot:P05067] |

**Important correction for feature-ID matching:** the study's feature IDs use the pattern `sp|<accession>|<ENTRY_NAME>`. For APP this is **`sp|P05067|A4_HUMAN`**, not `sp|P05067|APP_HUMAN` — UniProt's historical entry name for this protein is `A4_HUMAN` (from its original name "amyloid A4 protein"), and this has not changed despite the gene symbol being `APP`. Any lookup/parsing logic that assumes entry-name == gene-symbol-plus-`_HUMAN` will silently miss this protein. This was directly confirmed by fetching the live UniProt flat file for P05067 [uniprot:P05067].

**Biomarker forms and their position in the APP sequence** (UniProt `Contains:` sub-features, canonical 770-aa numbering) [uniprot:P05067]:

- **Amyloid-beta protein 42 (Abeta42)**: chain residues **672-713** (feature id `PRO_0000000092`).
- **Amyloid-beta protein 40 (Abeta40)**: chain residues **672-711** (feature id `PRO_0000000093`).

Both are proteolytic products of sequential beta- and gamma-secretase cleavage of APP; in a bottom-up (tryptic) MS search space they are not separate database entries but specific C-terminal fragments of the `A4_HUMAN`/APP sequence, distinguishable from each other and from the rest of APP only by the exact residue range and cleavage-site peptide(s) recovered — i.e., protein-level (PSM-to-`P05067`) identification of APP does **not** imply Abeta42 or Abeta40 was specifically observed.

**Plasma biomarker evidence, explicitly labeled by context:**

- **Abeta42/Abeta40 ratio and APP669-711/Abeta1-42 ratio — PLASMA/BLOOD, strong when measured by IP-MS.** Nakamura et al. (2018, *Nature*) measured these ratios in plasma from Japanese (discovery, n=121) and Australian (validation, n=252) cohorts and found high concordance with brain amyloid-PET status (AUC 96.7% discovery / 94.1% validation), explicitly using "immunoprecipitation coupled with mass spectrometry" [pmid:29420472]. Schindler et al. (2019, *Neurology*) independently replicated this with a WashU/C2N high-precision **immunoprecipitation + LC-MS/MS** assay (AUC 0.88 alone, 0.94 with age+APOE4) [pmid:31371569].
- **Analytically hard / small fold-change — PLASMA/BLOOD.** Weber et al. (2024, *Front Neurol*) report that even a purpose-optimized LC-MS/MS plasma Abeta42/40 assay achieved only an ~18% mean-ratio separation between amyloid-PET-positive and -negative groups — explicitly called larger than the ~8-11% separation typical of other published LC-MS/MS plasma Abeta42/40 assays — and note that 10% assay imprecision has been shown to produce a 26% misclassification rate for immunoassay-based Abeta42/40 [pmid:38595851]. This is the honest caveat the task asked for: plasma Abeta42/40 is a real, PET-concordant signal, but it is a **small-effect-size, high-precision-required** biomarker, unlike p-tau217/181, which show large fold changes.
- **PLASMA-EV specifically**: no plasma-EV-fraction-specific clinical-performance data for Abeta42/40 was located in this search. The only EV-mechanistic evidence found is cell-culture (HeLa/N2a) plus human AD brain-plaque tissue showing Abeta can be released in association with exosomes generated from multivesicular bodies [pmid:16837572] — mechanistically plausible grounds for why APP-derived material could appear in an EV prep, but **not** evidence of a validated blood/plasma-EV Abeta biomarker. This is flagged as an evidence gap, not a claim.
- **CSF/brain**: not the focus of this bounded topic; CSF Abeta42 and CSF Abeta42/40 are well-established (not separately verified here since out of scope).

### 2. Tau axis — parent protein MAPT

**Identity (verified against live UniProt record, 2026-09-15):**

| Field | Value | Source |
|---|---|---|
| Gene symbol | `MAPT` | UniProt GN line, ECO:0000312\|HGNC:HGNC:6893 [uniprot:P10636] |
| UniProt accession (human) | `P10636` | [uniprot:P10636] |
| UniProt entry name | **`TAU_HUMAN`** | [uniprot:P10636] |
| Canonical (displayed) isoform | "PNS-tau" (isoform 1, `P10636-1`), 758 aa | [uniprot:P10636] |
| Recommended full name | Microtubule-associated protein tau | [uniprot:P10636] |

This matches the study's expected feature ID exactly: **`sp|P10636|TAU_HUMAN`** [uniprot:P10636].

**Isoform-numbering caveat (important for any positional/phospho-site annotation):** UniProt's canonical, *displayed* sequence for P10636 is **PNS-tau**, a 758-residue isoform expressed in the peripheral nervous system. The AD biomarker literature (p-tau181, p-tau217, p-tau231) conventionally numbers phospho-sites against the **441-residue "Tau-F"/2N4R isoform** — the longest CNS isoform, which UniProt lists separately as isoform `P10636-8`, *not* as the canonical sequence [uniprot:P10636]. These two numbering systems are **not** interchangeable by simple offset arithmetic, because PNS-tau contains large inserted segments relative to Tau-441 (MAPT has 9 annotated isoforms differing by inclusion/exclusion of up to 5 of 15 exons) [uniprot:P10636]. To resolve this precisely, the canonical (`P10636`) and Tau-F (`P10636-8`) FASTA sequences were fetched directly from UniProt and aligned computationally: Tau-441 **Thr181 -> canonical Thr498**, Tau-441 **Thr217 -> canonical Thr534**, and Tau-441 **Thr231 -> canonical Thr548**. All three canonical positions (498, 534, 548) are independently annotated as `MOD_RES` **phosphothreonine** sites in the live UniProt feature table for P10636, which cross-validates the alignment-derived mapping [uniprot:P10636]. **Practical implication for this project:** if peptide-level phospho-site positions are ever reported against the UniProt canonical (758-aa) numbering rather than the literature's 441-aa numbering, p-tau181/217/231 will appear at residues 498/534/548, not 181/217/231 — a common source of silent mis-annotation.

**Plasma biomarker evidence per tau species, explicitly labeled by context:**

| Marker | Evidence context | Assay technology | Key plasma performance | Reference |
|---|---|---|---|---|
| Total tau (unmodified MAPT) | **PLASMA/BLOOD** — weak | Ultrasensitive immunoassay | Associated with AD at group level but weak, cohort-inconsistent correlations; large AD/normal-aging overlap; authors: "do not support plasma tau as an AD biomarker in individual people" | [pmid:27694257] |
| p-tau181 (Thr181) | **PLASMA/BLOOD** — strong, well-replicated | Ultrasensitive immunoassay (Simoa) | AUC 0.90-0.98 vs. cognitively unimpaired across 4 cohorts (n=1131); distinguishes AD from FTD/VaD/PSP-CBS/PD-MSA | [pmid:32333900] |
| p-tau217 (Thr217) | **PLASMA/BLOOD** — strongest single ATN plasma marker | Immunoassay | AUC 0.89-0.96 across 3 cohorts (n=1402); elevated ~20 years before estimated MCI onset in autosomal-dominant AD carriers; not significantly different from CSF/PET | [pmid:32722745] |
| p-tau231 (Thr231) | **PLASMA/BLOOD** — promising, earlier-rising, less replicated | Ultrasensitive immunoassay (Simoa) | AUC 0.92-0.94 vs. amyloid-negative CU; rises before p-tau181 and before amyloid-PET positivity threshold | [pmid:33585983] |
| p-tau217 & p-tau181 (MS-based) | **PLASMA/BLOOD**, but only after immunopurification | Immunoprecipitation + mass spectrometry | AUROC 0.99/0.98 (discovery, n=36), 0.92/0.75 (validation, n=92) — requires "immunopurified against N-terminus and mid-domain tau antibodies" before MS | [pmid:32725127] |

No CSF-only or plasma-EV-fraction-specific clinical-performance citation was sought or is claimed here for the tau markers (out of the bounded scope, and none of the assays above are EV-fraction-specific — all are whole/bulk plasma).

### 3. Practical note for plasma DIA-MS / the QC dynamic-range plot

This is the operationally important synthesis for annotating a QC dynamic-range plot from an **untargeted, plasma-EV DIA** experiment:

1. **The plasma proteome's dynamic range is the binding constraint.** More than 10 orders of magnitude separate albumin from the rarest clinically measured plasma analytes [pmid:12488461], and single-run, non-depleted/non-enriched LC-MS plasma profiling reproducibly reaches only a few hundred proteins (up to ~1,000 with extensive pre-fractionation), dominated by moderate-to-high-abundance, FDA-biomarker-class proteins [pmid:27135364] — not brain-/neuron-derived low-abundance proteins like APP-derived Abeta or MAPT-derived phospho-tau peptides.
2. **Every plasma study that successfully quantified these ATN markers by mass spectrometry used targeted immunoaffinity enrichment first** — immunoprecipitation coupled to MS (IP-MS) for Abeta42/40 [pmid:29420472; pmid:31371569] and for p-tau181/217 [pmid:32725127] — never untargeted/global shotgun or DIA proteomics on unenriched plasma. The remaining, most clinically validated plasma assays (p-tau181, p-tau217, p-tau231) are **ultrasensitive immunoassays (Simoa/equivalent)**, not mass spectrometry at all [pmid:32333900; pmid:32722745; pmid:33585983].
3. **Consequently, untargeted plasma-EV DIA-MS should not be expected to detect or quantify the diagnostic Abeta42/Abeta40 peptides or the p-tau181/217/231 phosphopeptides** at native plasma/EV abundance without a dedicated enrichment step (immunoaffinity, phospho-enrichment, or similar) — this is the correct, literature-grounded caveat to attach to a QC dynamic-range figure showing these features near or below the assay's practical detection floor.
4. **A protein-level ID is not a biomarker-level ID.** Because Abeta42/40 and p-tau181/217/231 are proteolytic fragments/phospho-forms of the single parent sequences `A4_HUMAN` (P05067) and `TAU_HUMAN` (P10636), a PSM or protein-group hit to `sp|P05067|A4_HUMAN` or `sp|P10636|TAU_HUMAN` in a DIA search only establishes that *some* peptide of the parent protein was observed — it does **not** by itself establish that the Abeta42/40 C-terminal fragment or a Thr181/217/231-phosphorylated peptide specifically was observed. That must be checked at the individual peptide/phosphopeptide level (and, for MAPT, against the correct isoform-numbering convention noted above) before any ATN-marker claim is attached to a feature ID in this dataset.
5. **If APP or MAPT peptides are nonetheless observed in the plasma-EV DIA data**, the most literature-plausible explanation is generic/structural tryptic peptides from the bulk parent protein (APP is a type-I transmembrane protein and a plausible generic EV/membrane-cargo protein; cf. the exosome-association mechanism reported in cell culture and brain tissue [pmid:16837572]) rather than the specific low-abundance diagnostic Abeta or phospho-tau epitopes, absent enrichment.

## Methods / sources consulted

- **UniProt** (`rest.uniprot.org`): fetched the live flat-file records (`.txt`) and canonical/isoform FASTA sequences for `P05067` and `P10636` directly (2026-09-15) to confirm accession, entry name, gene symbol, organism, canonical sequence length, isoform list, and `MOD_RES`/`CHAIN` feature annotations. Performed an independent computational sequence alignment (Python) between the UniProt canonical P10636 sequence and isoform `P10636-8` (Tau-F/2N4R, 441 aa) to resolve the Thr181/217/231 literature-numbering-to-canonical-numbering mapping, then cross-checked the resulting canonical positions (498, 534, 548) against UniProt's own `MOD_RES` phosphothreonine annotations at those exact positions — an independent, self-consistent verification.
- **PubMed / NCBI E-utilities** (`eutils.ncbi.nlm.nih.gov`): used `esearch`/`esummary`/`efetch` to identify and pull full abstracts for all cited primary literature, confirming title, journal, year, volume/issue/pages, DOI, and PMID for every reference before citing it (i.e., every PMID/DOI below was independently resolved and its abstract read, not assumed from a web-search snippet).
- **CrossRef API** (`api.crossref.org`): used to cross-check DOIs for two references (Mattsson et al. 2016; Anderson & Anderson 2002) where the PubMed summary record's `elocationid` was ambiguous or absent; the PubMed-`efetch`-confirmed DOI was used as authoritative in the final reference where the two differed.
- **Web search**: used only for discovery/triage (identifying candidate papers and general context), never as the basis for a citation — every citation used here was subsequently verified directly against PubMed (title/journal/year/DOI match) or, for the Weber et al. 2024 fold-change numbers, against the article's PMC full text via a targeted fetch to obtain and quote the exact sentences.
- **PMC full text** (`pmc.ncbi.nlm.nih.gov`): fetched directly to extract and quote the exact fold-change/misclassification-rate figures for the plasma Abeta42/40 LC-MS/MS analytical-precision claim (Weber et al. 2024), since these numbers were not present in the PubMed abstract.

**Claims not independently verifiable / explicitly flagged as gaps:**
- No plasma-**EV-fraction-specific** (as opposed to whole/bulk plasma) clinical-performance literature for Abeta42/40 or p-tau181/217/231 was located; this is stated as an evidence gap in the Summary and Detailed Findings, not asserted as fact.
- The specific analytical platform used in Palmqvist et al. 2020 (JAMA) for plasma P-tau217 (e.g., Meso Scale Discovery electrochemiluminescence vs. another electrochemiluminescence/immunoassay format) is not stated in the verified abstract and was deliberately left as generic "immunoassay" rather than guessed.
- Typical absolute plasma concentrations (pg/mL) for Abeta42, Abeta40, and the p-tau species were not included as numeric claims in this document because a specific, independently verifiable primary-source number for each was not confirmed during this research pass; the qualitative "far below the untargeted plasma-proteome detection floor, requires ultrasensitive immunoassay or IP-MS" conclusion is instead supported by the assay-technology evidence in references [pmid:29420472], [pmid:31371569], [pmid:32725127], [pmid:32333900], [pmid:32722745], [pmid:33585983] taken together with the general plasma dynamic-range literature [pmid:12488461; pmid:27135364].

## References

1. UniProt Consortium. UniProtKB entry **P05067** (APP_HUMAN record; entry name `A4_HUMAN`), Homo sapiens Amyloid-beta precursor protein. https://www.uniprot.org/uniprotkb/P05067/entry (accessed 2026-09-15, live REST record).
2. UniProt Consortium. UniProtKB entry **P10636** (`TAU_HUMAN`), Homo sapiens Microtubule-associated protein tau. https://www.uniprot.org/uniprotkb/P10636/entry (accessed 2026-09-15, live REST record; isoform P10636-8 "Tau-F"/2N4R also consulted).
3. Nakamura A, Kaneko N, Villemagne VL, et al. High performance plasma amyloid-β biomarkers for Alzheimer's disease. *Nature*. 2018;554(7691):249-254. doi:10.1038/nature25456. PMID:29420472.
4. Schindler SE, Bollinger JG, Ovod V, et al. High-precision plasma β-amyloid 42/40 predicts current and future brain amyloidosis. *Neurology*. 2019;93(17):e1647-e1659. doi:10.1212/WNL.0000000000008081. PMID:31371569.
5. Weber DM, Taylor SW, Lagier RJ, et al. Clinical utility of plasma Aβ42/40 ratio by LC-MS/MS in Alzheimer's disease assessment. *Front Neurol*. 2024;15:1364658. doi:10.3389/fneur.2024.1364658. PMID:38595851.
6. Rajendran L, Honsho M, Zahn TR, et al. Alzheimer's disease beta-amyloid peptides are released in association with exosomes. *Proc Natl Acad Sci U S A*. 2006;103(30):11172-11177. doi:10.1073/pnas.0603838103. PMID:16837572.
7. Karikari TK, Pascoal TA, Ashton NJ, et al. Blood phosphorylated tau 181 as a biomarker for Alzheimer's disease: a diagnostic performance and prediction modelling study using data from four prospective cohorts. *Lancet Neurol*. 2020;19(5):422-433. doi:10.1016/S1474-4422(20)30071-5. PMID:32333900.
8. Palmqvist S, Janelidze S, Quiroz YT, et al. Discriminative Accuracy of Plasma Phospho-tau217 for Alzheimer Disease vs Other Neurodegenerative Disorders. *JAMA*. 2020;324(8):772-781. doi:10.1001/jama.2020.12134. PMID:32722745.
9. Ashton NJ, Pascoal TA, Karikari TK, et al. Plasma p-tau231: a new biomarker for incipient Alzheimer's disease pathology. *Acta Neuropathol*. 2021;141(5):709-724. doi:10.1007/s00401-021-02275-6. PMID:33585983.
10. Mattsson N, Zetterberg H, Janelidze S, et al. Plasma tau in Alzheimer disease. *Neurology*. 2016;87(17):1827-1835. doi:10.1212/WNL.0000000000003246. PMID:27694257.
11. Barthélemy NR, Horie K, Sato C, Bateman RJ. Blood plasma phosphorylated-tau isoforms track CNS change in Alzheimer's disease. *J Exp Med*. 2020;217(11):e20200861. doi:10.1084/jem.20200861. PMID:32725127.
12. Anderson NL, Anderson NG. The human plasma proteome: history, character, and diagnostic prospects. *Mol Cell Proteomics*. 2002;1(11):845-867. doi:10.1074/mcp.r200007-mcp200. PMID:12488461.
13. Geyer PE, Kulak NA, Pichler G, Holdt LM, Teupser D, Mann M. Plasma Proteome Profiling to Assess Human Health and Disease. *Cell Syst*. 2016;2(3):185-195. doi:10.1016/j.cels.2016.02.015. PMID:27135364.
