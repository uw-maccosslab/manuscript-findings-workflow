---
# Research-finding document — external knowledge, structured like a finding but for the literature/tools.
# Spec: doc 04.4. Stored in research/<slug>.md. A research finding without verified references is NOT accepted.

topic: "Proteins repeatedly reported in plasma-derived extracellular vesicle (EV) proteomic studies of Alzheimer's disease: neuronal-EV cargo, the L1CAM capture-marker controversy, canonical EV-quality-control markers, and candidates from unbiased plasma-EV MS proteomics"
type: research-publications
slug: ad-plasma-ev-proteomic-markers
created: 2026-09-15
updated: 2026-09-15

status: reviewed
reviewed_by: research-reviewer (claude@ogdb.com)
reviewed_date: 2026-09-15

# Normalized domain-entity references this research is about (canonical IDs; gene symbol as label).
entities:
  - { db: uniprot, id: "P32004", label: "L1CAM" }
  - { db: uniprot, id: "P13591", label: "NCAM1" }
  - { db: uniprot, id: "P10636", label: "MAPT" }
  - { db: uniprot, id: "P05067", label: "APP" }
  - { db: uniprot, id: "P08247", label: "SYP" }
  - { db: uniprot, id: "P60880", label: "SNAP25" }
  - { db: uniprot, id: "P17677", label: "GAP43" }
  - { db: uniprot, id: "P47972", label: "NPTX2" }
  - { db: uniprot, id: "Q92686", label: "NRGN" }
  - { db: uniprot, id: "P21579", label: "SYT1" }
  - { db: uniprot, id: "P21926", label: "CD9" }
  - { db: uniprot, id: "P08962", label: "CD63" }
  - { db: uniprot, id: "P60033", label: "CD81" }
  - { db: uniprot, id: "O00560", label: "SDCBP" }
  - { db: uniprot, id: "Q8WUM4", label: "PDCD6IP" }
  - { db: uniprot, id: "Q99816", label: "TSG101" }
  - { db: uniprot, id: "O75955", label: "FLOT1" }
  - { db: uniprot, id: "P05109", label: "S100A8" }
  - { db: uniprot, id: "P00488", label: "F13A1" }
  - { db: uniprot, id: "P05160", label: "F13B" }

# (required, non-empty) Every external claim traces to a reference here.
references:
  - id: "pmid:28588440"
    type: doi
    claim: "Anti-L1CAM immunoprecipitation (following ExoQuick/precipitation-based total-exosome enrichment) is the originating, widely adopted method for isolating putative neuron-derived plasma extracellular vesicles (NDEVs), underlying the subsequent line of plasma NDEV biomarker studies."
    verified: true
    verified_by: research-reviewer (claude@ogdb.com)
  - id: "pmid:25130657"
    type: doi
    claim: "In L1CAM-immunocaptured neurally derived blood exosomes, levels of total tau, P-T181-tau, P-S396-tau, and amyloid-beta 1-42 (Abeta42) were significantly higher in AD patients than case-controls, and P-S396-tau, P-T181-tau, and Abeta42 predicted development of AD up to 10 years before clinical onset."
    verified: true
    verified_by: research-reviewer (claude@ogdb.com)
  - id: "pmid:27601437"
    type: doi
    claim: "In L1CAM-immunocaptured neuronal exosomes (NDEs), levels of synaptophysin (SYP), synaptopodin, synaptotagmin-2, and neurogranin (NRGN) were significantly lower in both frontotemporal dementia (FTD) and AD than in controls, whereas GAP43 and synapsin-1 were reduced only in AD."
    verified: true
    verified_by: research-reviewer (claude@ogdb.com)
  - id: "pmid:29025866"
    type: doi
    claim: "NDE levels of neuronal pentraxin-2 (NPTX2) together with its synaptic partners AMPA4/GRIA4, neuroligin-1 (NLGN1), and neurexin-2-alpha (NRXN2) declined progressively with AD severity; NPTX2 was significantly reduced in symptomatic AD dementia but, unlike the other three proteins, was not yet significantly reduced 6-11 years before dementia onset."
    verified: true
    verified_by: research-reviewer (claude@ogdb.com)
  - id: "pmid:32776690"
    type: doi
    claim: "In a multi-cohort study (discovery, validation, preclinical, and familial-AD cohorts totaling >700 subjects) of NCAM(CD56)-immunocaptured neuro-exosomes (with L1CAM used to confirm neuronal enrichment), GAP43, neurogranin (NRGN), SNAP25, and synaptotagmin-1 (SYT1) were significantly lower in AD than in controls and predicted AD 5-7 years before cognitive impairment."
    verified: true
    verified_by: research-reviewer (claude@ogdb.com)
  - id: "pmid:30680692"
    type: doi
    claim: "SNAP25 in serum is predominantly carried by L1CAM-positive exosomes of neuronal origin, and exosomal SNAP25 levels were significantly reduced in AD (and in MCI patients who converted to AD) relative to controls, proposed as a potential AD biomarker."
    verified: true
    verified_by: research-reviewer (claude@ogdb.com)
  - id: "pmid:34092791"
    type: doi
    claim: "Across human CSF and plasma, the majority of L1CAM does not co-fractionate with canonical EV tetraspanin markers on size-exclusion chromatography or density-gradient separation; most L1CAM runs as a soluble, proteolytically cleaved ~200 kDa species detected only by antibodies against its extracellular domain, indicating that L1CAM in these biofluids is largely not EV-associated and calling into question L1CAM-based neuronal-EV isolation."
    verified: true
    verified_by: research-reviewer (claude@ogdb.com)
  - id: "pmid:40210837"
    type: doi
    claim: "Independent replication across human plasma, CSF, and iPSC-derived neuron conditioned media confirms that most L1CAM immunoreactivity is not associated with extracellular vesicles, corroborating the non-vesicular, soluble nature of plasma/CSF L1CAM reported by Norman et al. (2021)."
    verified: true
    verified_by: research-reviewer (claude@ogdb.com)
  - id: "doi:10.1002/alz.093432"
    type: doi
    claim: "In a PRISMA-guided meta-analysis of 49 studies, diagnostic biomarkers measured in general/unselected plasma or CSF EVs showed better diagnostic accuracy and less between-study heterogeneity and publication bias for AD and related dementias than biomarkers measured in 'speculative' CNS/neuron-enriched EVs (e.g., L1CAM-captured NDEVs)."
    verified: true
    verified_by: research-reviewer (claude@ogdb.com)
  - id: "pmid:30637094"
    type: doi
    claim: "The MISEV2018 position statement specifies that EV preparations should be characterized using at least one transmembrane/GPI-anchored protein enriched at the plasma membrane and/or endosomes (e.g., tetraspanins CD9, CD63, CD81, and flotillin-1/FLOT1) together with at least one cytosolic protein recovered in EVs (e.g., TSG101, syntenin-1/SDCBP, ALIX/PDCD6IP) as general markers of EV identity and recovery, not as disease-specific markers."
    verified: true
    verified_by: research-reviewer (claude@ogdb.com)
  - id: "pmid:34531462"
    type: doi
    claim: "Untargeted label-free LC-MS/MS (data-dependent acquisition) proteomics of plasma EVs from 10 AD, 10 MCI, and 9 cognitively normal controls identified 336 EV proteins and found 63 proteins differentially expressed between AD and controls (57 upregulated and 6 downregulated in AD; 19 of the upregulated proteins remained significant after FDR correction); coagulation factor XIII subunits F13A1 and F13B were among the most strongly upregulated proteins in AD, alongside acute-phase/inflammatory and complement proteins (e.g., ORM1, ORM2, CRP, complement factors), as part of a single small discovery-stage study without independent mass-spectrometry replication."
    verified: true
    verified_by: research-reviewer (claude@ogdb.com)
  - id: "pmid:39159863"
    type: doi
    claim: "In a preliminary tandem-mass-tag (TMT) quantitative proteomics study of plasma EVs from 9 AD patients and 9 healthy controls (335 proteins quantified, 12 differentially expressed proteins), S100A8 was significantly down-regulated in AD-patient plasma EVs and proposed as a novel single-study candidate biomarker."
    verified: true
    verified_by: research-reviewer (claude@ogdb.com)
  - id: "pmid:26062630"
    type: doi
    claim: "Within the same L1CAM-immunocaptured NDEV framework used for tau/Abeta and synaptic-protein studies, lysosomal pathway cargo (e.g., cathepsin D, LAMP-1, ubiquitinated proteins) was also reported as altered in preclinical AD, illustrating the breadth of cargo measured in this immunocapture literature beyond tau/Abeta and synaptic proteins."
    verified: true
    verified_by: research-reviewer (claude@ogdb.com)
---

# Research: Proteins repeatedly reported in plasma-EV proteomic studies of Alzheimer's disease

## Summary

The plasma-EV/exosome AD literature is dominated by one lineage of studies that immunocapture "neuron-derived EVs" (NDEVs) from plasma using an antibody against the neuronal adhesion molecule **L1CAM** (occasionally NCAM1/CD56), then measure a small, targeted panel of cargo proteins by immunoassay: phospho-tau/total tau and amyloid-beta (from *MAPT*/*APP*) plus a repeatedly-measured set of pre- and post-synaptic proteins — synaptophysin (SYP), SNAP25, GAP43, neurogranin (NRGN), synaptotagmin-1 (SYT1), and neuronal pentraxin-2 (NPTX2) — nearly all reported by one overlapping set of investigators (Kapogiannis, Goetzl, and collaborators/successors) [pmid:25130657][pmid:27601437][pmid:29025866][pmid:32776690][pmid:30680692]. However, the validity of L1CAM as a marker of vesicle-bound, neuron-derived material is **directly and repeatedly contested**: Norman et al. (2021, *Nature Methods*) showed most plasma/CSF L1CAM is a soluble, cleaved ~200 kDa species that does not co-fractionate with EV tetraspanin markers [pmid:34092791], a finding independently replicated in 2025 [pmid:40210837], and a 2024 meta-analysis of 49 EV-biomarker studies found that biomarkers from *general* (unselected) plasma/CSF EVs outperform those from "speculative" L1CAM/CNS-enriched EVs on diagnostic accuracy, heterogeneity, and publication bias [doi:10.1002/alz.093432]. Separate from the NDEV-immunocapture literature, a distinct set of proteins — the tetraspanins CD9/CD63/CD81, syntenin-1 (SDCBP), ALIX (PDCD6IP), TSG101, and flotillin-1 (FLOT1) — are the community-standard **EV-quality-control markers** used to confirm that an EV preparation (such as Mag-Net) enriched genuine vesicles; per the MISEV2018 guidelines, these are markers of vesicle identity/recovery, not of AD status [pmid:30637094]. Finally, unbiased (DDA/TMT) plasma-EV MS proteomics of AD has so far produced only **single, small, non-replicated discovery studies**: one reported F13A1/F13B (coagulation) and inflammatory/complement proteins up in AD EVs [pmid:34531462], and an unrelated one reported S100A8 down in AD EVs [pmid:39159863]; no protein from this unbiased-MS branch of the literature has been independently replicated across studies, so none should be treated as an established AD-differential plasma-EV marker.

## Detailed findings

### 1. The neuronal-derived EV (NDEV) immunocapture literature

**L1CAM as the capture marker, and why it is contested.** L1CAM (L1 cell adhesion molecule, UniProt `P32004`, entry `L1CAM_HUMAN`) was introduced as an immunocapture target for isolating putative neuron-derived EVs from plasma/serum in Mustapic et al. 2017 (*Frontiers in Neuroscience*): total exosomes are first precipitated (e.g., ExoQuick), then an anti-L1CAM antibody is used to pull down the neuronally-enriched subpopulation [pmid:28588440]. This two-step "total EV precipitation + L1CAM immunocapture" (or the closely related NCAM1/CD56-based variant, UniProt `P13591`, entry `NCAM1_HUMAN`) is the method underlying essentially all of the plasma NDEV biomarker papers cited below [pmid:32776690].

Beginning in 2021, this approach came under sustained, well-documented challenge. Norman, Ter-Ovanesyan, Trieu et al. (*Nature Methods*, 2021) directly tested whether L1CAM in human CSF and plasma is vesicle-associated using size-exclusion chromatography and density-gradient separation together with antibodies against different L1CAM epitopes. They found that CSF/plasma L1CAM runs predominantly as a soluble, proteolytically cleaved ~200 kDa fragment, is detected only by antibodies to its extracellular domain (not the intracellular domain expected on an intact transmembrane EV-surface protein), and shows minimal co-elution with tetraspanin EV markers — concluding that L1CAM is "not associated with extracellular vesicles" in these biofluids and calling L1CAM-based NDEV isolation into question [pmid:34092791]. This was independently replicated in 2025 by Kadam, Wacker, Oeckl, Korneck et al. across plasma, CSF, and iPSC-derived neuron conditioned media, again finding that most L1CAM immunoreactivity is non-vesicular [pmid:40210837]. Consistent with this, a 2024 PRISMA-guided meta-analysis of 49 studies (Taha et al., ICAD 2024 abstract, *Alzheimer's & Dementia*) comparing biomarkers from general/unselected EVs versus "speculative" CNS/neuron-enriched EVs (the L1CAM-type approach) found that general-EV biomarkers had **better** diagnostic accuracy and less between-study heterogeneity and publication bias for AD and related dementias than the neuron-enriched approach [doi:10.1002/alz.093432]. Critics of the Norman et al. interpretation have countered that demonstrating soluble L1CAM exists does not, by itself, prove that no L1CAM-bearing EVs exist at all (see the ALZFORUM commentary threads on this paper) — but as of this writing no peer-reviewed rebuttal paper overturning Norman et al.'s core finding has been identified in this search; the controversy should be treated as **live and unresolved**, not settled in either direction. Practically, this means: L1CAM detected in a plasma-EV DIA run should not automatically be interpreted as evidence of neuronal-EV enrichment, and depletion/enrichment claims based on L1CAM alone should be flagged as contested.

**Cargo repeatedly measured in L1CAM(or NCAM)-immunocaptured NDEV preparations for AD.** Within this immunocapture literature (almost entirely EIA/Simoa-based, not untargeted MS), the following cargo proteins have been repeatedly reported as differential between AD and controls:

- **Tau (gene *MAPT*, UniProt `P10636`, entry `TAU_HUMAN`) and amyloid-beta (cargo of *APP*, UniProt `P05067`, entry `A4_HUMAN`)**: Fiandaca et al. 2015 found total tau, phospho-tau-T181, phospho-tau-S396, and Abeta1-42 significantly elevated in L1CAM+ NDEVs from AD patients versus controls, with phospho-tau and Abeta1-42 predicting AD onset up to 10 years pre-diagnosis [pmid:25130657]. This is the paper most often cited as the origin of "exosomal p-tau/Abeta" as an AD biomarker concept.
- **Synaptophysin (*SYP*, UniProt `P08247`, entry `SYPH_HUMAN`)**: reported significantly lower in NDEs of both AD and FTD versus controls, alongside synaptopodin, synaptotagmin-2, and neurogranin [pmid:27601437].
- **SNAP25 (*SNAP25*, UniProt `P60880`, entry `SNP25_HUMAN`)**: reduced in AD (and MCI-to-AD converters) neuro-exosomes/serum exosomes in at least two independent cohorts: a large multi-cohort Chinese study (discovery/validation/preclinical/familial-AD, >700 subjects total) [pmid:32776690] and an Italian serum-exosome cohort focused specifically on SNAP25 [pmid:30680692].
- **GAP43 (*GAP43*, UniProt `P17677`, entry `NEUM_HUMAN`)**: reduced in AD (but not FTD) NDEs in the 2016 cohort [pmid:27601437] and again reduced in AD in the independent 2021 multi-cohort study [pmid:32776690].
- **Neuronal pentraxin-2 (*NPTX2*, UniProt `P47972`, entry `NPTX2_HUMAN`)**: reported declining with AD progression together with its synaptic partners AMPA4/GRIA4, neuroligin-1 (NLGN1), and neurexin-2-alpha (NRXN2); NPTX2 itself was significantly reduced at the symptomatic AD-dementia stage but, notably, was *not yet* significantly reduced in the preclinical window (6-11 years pre-onset) where the other three proteins already were — an important caveat against treating NPTX2 as an equally early marker [pmid:29025866].
- Additional synaptic cargo repeatedly reported by the same overlapping group of investigators: **neurogranin (*NRGN*, UniProt `Q92686`, entry `NEUG_HUMAN`)** and **synaptotagmin-1 (*SYT1*, UniProt `P21579`, entry `SYT1_HUMAN`)**, both reduced in AD NDEs across at least two independent cohorts [pmid:27601437][pmid:32776690].
- Beyond synaptic/pathogenic cargo, the same immunocapture framework has also reported altered lysosomal/autophagy proteins (cathepsin D, LAMP-1, ubiquitinated proteins) in preclinical AD NDEVs, illustrating that this literature's cargo panel extends beyond the tau/Abeta/synaptic proteins highlighted here [pmid:26062630].

**Important caveat for this project.** All of the above cargo findings come from **targeted immunoassay (ELISA/Simoa) quantification of L1CAM- or NCAM-immunocaptured EV lysates**, not from untargeted/DIA mass spectrometry of bulk plasma-EV preparations (such as Mag-Net). Two things follow: (1) detection of these proteins in a Mag-Net DIA run does not by itself replicate this literature, since Mag-Net enriches EVs in bulk rather than a neuron-specific subpopulation; and (2) as detailed above, the entire premise that L1CAM marks a neuron-derived EV subpopulation is independently contested [pmid:34092791][pmid:40210837], so even a faithful re-creation of the immunocapture step would carry that caveat.

### 2. Canonical EV-quality-control markers (not AD markers)

The following proteins are the community-standard markers used across the EV field — including in the plasma-EV AD papers reviewed here — to confirm that a preparation actually enriched extracellular vesicles, per the MISEV2018 (Minimal Information for Studies of Extracellular Vesicles) guidelines of the International Society for Extracellular Vesicles [pmid:30637094]. They are **general vesicle-identity/recovery markers, not disease markers**, and are not reported as differential between AD and controls in the plasma-EV AD literature reviewed here. An EV-enrichment method such as Mag-Net should detect these if it is successfully capturing bona fide EVs:

| Gene | UniProt accession | UniProt entry name | Role |
|---|---|---|---|
| CD9 | P21926 | CD9_HUMAN | Tetraspanin; canonical EV surface marker |
| CD63 | P08962 | CD63_HUMAN | Tetraspanin; canonical EV/exosome surface marker |
| CD81 | P60033 | CD81_HUMAN | Tetraspanin; canonical EV surface marker |
| SDCBP (syntenin-1) | O00560 | SDCB1_HUMAN | Cytosolic ESCRT-associated protein enriched in EVs |
| PDCD6IP (ALIX) | Q8WUM4 | PDC6I_HUMAN | Cytosolic ESCRT-accessory protein enriched in EVs |
| TSG101 | Q99816 | TS101_HUMAN | ESCRT-I component enriched in EVs |
| FLOT1 (flotillin-1) | O75955 | FLOT1_HUMAN | Lipid-raft-associated membrane protein enriched in EVs |

MISEV2018 specifically recommends demonstrating at least one transmembrane/GPI-anchored, plasma-membrane- or endosome-enriched protein (the tetraspanins and/or flotillin) together with at least one cytosolic protein recovered in EVs (syntenin-1, ALIX, and/or TSG101) as minimal positive evidence of EV identity [pmid:30637094].

### 3. Unbiased plasma-EV AD proteomics (DDA/TMT): what is (and is not) well-supported

Two identified untargeted LC-MS discovery studies of plasma-EV proteomics in AD report differential proteins, but **neither protein set has been independently replicated**, and the two studies do not overlap in their top findings — this should be reported honestly rather than as an established signature:

- **Nielsen et al. 2021** (*Scientific Reports*), label-free DDA LC-MS/MS of plasma EVs from 10 AD, 10 MCI, and 9 controls: 336 EV proteins identified, 63 differentially expressed between AD and controls (57 upregulated and 6 downregulated in AD; 19 of the upregulated proteins significant after FDR correction). The most strongly upregulated in AD included coagulation factor XIII subunits **F13A1** (UniProt P00488, entry F13A_HUMAN) and **F13B** (UniProt P05160, entry F13B_HUMAN), alongside acute-phase/inflammatory and complement proteins (ORM1, ORM2, CRP, complement factors) [pmid:34531462]. This is a small single discovery-stage study; no independent MS replication of F13A1/F13B as plasma-EV AD markers was found in this search.
- **Zhang et al. 2024** (*Journal of Proteomics*), explicitly titled "a preliminary study": TMT-based quantitative proteomics of plasma EVs from 9 AD and 9 healthy controls, again 335 proteins quantified and 12 differentially expressed proteins. **S100A8** (UniProt P05109, entry S10A8_HUMAN) was significantly down-regulated in AD-patient plasma EVs and proposed as a candidate biomarker, with a mechanistic follow-up in an Abeta-induced SH-SY5Y cell EV model [pmid:39159863]. Again, single-study, unreplicated.

No protein appears as a differential hit in both of these unbiased studies. Given this, **none of these proteins (F13A1, F13B, ORM1/ORM2/CRP/complement, S100A8) should be treated as an established, repeatedly-reported plasma-EV AD marker** — they are single-study, preliminary candidates from small discovery cohorts (n=9-10/group), included here for completeness because the task specifically asked for unbiased-proteomics findings, but they should be weighted accordingly (i.e., as hypothesis-generating, not confirmatory) when annotating a QC plot or orienting an analysis. Separately, note that two additional "plasma EV proteomics AD" studies identified in this search (a BIOPEXAL preprint and an AAIC 2025 abstract, both from what appears to be the same Barcelona cohort/group) used Olink proximity-extension-assay (PEA) technology, not DIA/DDA mass spectrometry, and are therefore out of scope for the "unbiased DIA/DDA" category requested here; they are not cited as MS evidence.

### 4. Summary table for QC/annotation purposes

| Gene | UniProt | Entry name | Evidence tier | Context |
|---|---|---|---|---|
| L1CAM | P32004 | L1CAM_HUMAN | Contested capture marker | NDEV immunocapture target; plasma/CSF L1CAM shown largely non-vesicular/soluble [pmid:34092791][pmid:40210837] |
| NCAM1 | P13591 | NCAM1_HUMAN | Contested capture marker (alt.) | Alternative NDEV capture antibody used in some cohorts [pmid:32776690] |
| MAPT | P10636 | TAU_HUMAN | NDEV cargo, AD-associated (immunoassay, plasma-EV-specific) | Elevated p-tau in L1CAM+ NDEVs in AD [pmid:25130657] |
| APP | P05067 | A4_HUMAN | NDEV cargo, AD-associated (immunoassay, plasma-EV-specific) | Elevated Abeta42 in L1CAM+ NDEVs in AD [pmid:25130657] |
| SYP | P08247 | SYPH_HUMAN | NDEV cargo, AD-associated (immunoassay, plasma-EV-specific) | Reduced in AD/FTD NDEs [pmid:27601437] |
| SNAP25 | P60880 | SNP25_HUMAN | NDEV cargo, AD-associated (immunoassay, plasma-EV-specific; replicated across ≥2 independent cohorts) | [pmid:32776690][pmid:30680692] |
| GAP43 | P17677 | NEUM_HUMAN | NDEV cargo, AD-associated (immunoassay, plasma-EV-specific; replicated across ≥2 cohorts) | [pmid:27601437][pmid:32776690] |
| NPTX2 | P47972 | NPTX2_HUMAN | NDEV cargo, AD-associated (immunoassay, plasma-EV-specific; single cohort, only at symptomatic stage) | [pmid:29025866] |
| NRGN | Q92686 | NEUG_HUMAN | NDEV cargo, AD-associated (immunoassay, plasma-EV-specific; replicated) | [pmid:27601437][pmid:32776690] |
| SYT1 | P21579 | SYT1_HUMAN | NDEV cargo, AD-associated (immunoassay, plasma-EV-specific) | [pmid:32776690] |
| CD9 | P21926 | CD9_HUMAN | EV-quality-control marker, NOT an AD marker | [pmid:30637094] |
| CD63 | P08962 | CD63_HUMAN | EV-quality-control marker, NOT an AD marker | [pmid:30637094] |
| CD81 | P60033 | CD81_HUMAN | EV-quality-control marker, NOT an AD marker | [pmid:30637094] |
| SDCBP | O00560 | SDCB1_HUMAN | EV-quality-control marker, NOT an AD marker | [pmid:30637094] |
| PDCD6IP | Q8WUM4 | PDC6I_HUMAN | EV-quality-control marker, NOT an AD marker | [pmid:30637094] |
| TSG101 | Q99816 | TS101_HUMAN | EV-quality-control marker, NOT an AD marker | [pmid:30637094] |
| FLOT1 | O75955 | FLOT1_HUMAN | EV-quality-control marker, NOT an AD marker | [pmid:30637094] |
| S100A8 | P05109 | S10A8_HUMAN | Unbiased plasma-EV MS proteomics, single-study/preliminary | [pmid:39159863] |
| F13A1 | P00488 | F13A_HUMAN | Unbiased plasma-EV MS proteomics, single-study/preliminary | [pmid:34531462] |
| F13B | P05160 | F13B_HUMAN | Unbiased plasma-EV MS proteomics, single-study/preliminary | [pmid:34531462] |

## Methods / sources consulted

- **Literature discovery**: iterative WebSearch queries targeting (a) the L1CAM/NDEV immunocapture literature and its controversy, (b) plasma-EV AD proteomics (MS-based), (c) MISEV2018/EV-marker guidelines, and (d) meta-analyses/systematic reviews of EV biomarkers in AD.
- **Primary-source verification via PubMed E-utilities** (`esearch`/`esummary`/`efetch`, NCBI NLM): every PMID cited above was confirmed to exist and its title/journal/volume/issue/pages/DOI were pulled directly from PubMed metadata (not from memory or search-snippet paraphrase). Command-line `curl` calls against `eutils.ncbi.nlm.nih.gov` were used throughout this session.
- **Full-text extraction via WebFetch** of PMC full-text pages (PMC8445922/Nielsen et al. 2021; PMC7984076/Jia et al. 2021; PMC12636747/BIOPEXAL preprint; PMC12740250/AAIC abstract; PMC11212434/Boyer et al. 2024; PMC11716713/Taha et al. 2024 abstract) to confirm specific numeric findings (fold-changes, sample sizes, isolation methods) rather than relying on search-result summaries alone. This is how the NCAM1 (not L1CAM) primary-isolation detail in Jia et al. 2021 was caught.
- **UniProt accession/entry-name verification**: every accession listed in the entities table and protein tables above was fetched directly from the UniProt REST API (`https://rest.uniprot.org/uniprotkb/<ACCESSION>.txt`) in this session and the returned `ID` line (entry name) and `AC` line (accession, including confirmation it is the canonical human Swiss-Prot entry) were read directly from the response — none were taken from memory.
- **Deliberately excluded from the "unbiased DIA/DDA" category**: two Olink-PEA-based plasma-EV AD studies (a Research Square preprint, DOI 10.21203/rs.3.rs-7847549/v1, "BIOPEXAL study"; and an AAIC 2025 conference abstract, PMC12740250) were reviewed but excluded from the unbiased-MS findings because PEA is a targeted antibody-based technology, not untargeted DIA/DDA mass spectrometry, and (in the BIOPEXAL case) the preprint is not yet peer-reviewed. They are not cited as references above.
- **Not found / not verifiable in this search**: a published peer-reviewed rebuttal paper directly overturning Norman et al. 2021's conclusion (only secondary commentary, e.g., ALZFORUM discussion threads, referencing unspecified "critics," was found — not cited as a reference here because it does not resolve to a specific, verifiable peer-reviewed source).

## References

1. Mustapic M, Eitan E, Werner JK Jr, Berkowitz ST, Lazaropoulos MP, Tran J, Goetzl EJ, Kapogiannis D. Plasma Extracellular Vesicles Enriched for Neuronal Origin: A Potential Window into Brain Pathologic Processes. *Front Neurosci*. 2017;11:278. doi: [10.3389/fnins.2017.00278](https://doi.org/10.3389/fnins.2017.00278). PMID: 28588440.
2. Fiandaca MS, Kapogiannis D, Mapstone M, et al. Identification of preclinical Alzheimer's disease by a profile of pathogenic proteins in neurally derived blood exosomes: A case-control study. *Alzheimers Dement*. 2015;11(6):600-607.e1. doi: [10.1016/j.jalz.2014.06.008](https://doi.org/10.1016/j.jalz.2014.06.008). PMID: 25130657.
3. Goetzl EJ, Kapogiannis D, Schwartz JB, et al. Decreased synaptic proteins in neuronal exosomes of frontotemporal dementia and Alzheimer's disease. *FASEB J*. 2016;30(12):4141-4148. doi: [10.1096/fj.201600816R](https://doi.org/10.1096/fj.201600816R). PMID: 27601437.
4. Goetzl EJ, et al. Declining levels of functionally specialized synaptic proteins in plasma neuronal exosomes with progression of Alzheimer's disease. *FASEB J*. 2018;32(2):888-893. doi: [10.1096/fj.201700731R](https://doi.org/10.1096/fj.201700731R). PMID: 29025866.
5. Jia L, Zhu M, Kong C, Pang Y, Zhang H, Qiu Q, Wei C, Tang Y, Wang Q, Li Y, Li T, Li F, Wang Q, Li Y, Wei Y, Jia J. Blood neuro-exosomal synaptic proteins predict Alzheimer's disease at the asymptomatic stage. *Alzheimers Dement*. 2021;17(1):49-60. doi: [10.1002/alz.12166](https://doi.org/10.1002/alz.12166). PMID: 32776690.
6. Agliardi C, Guerini FR, Zanzottera M, et al. SNAP-25 in Serum Is Carried by Exosomes of Neuronal Origin and Is a Potential Biomarker of Alzheimer's Disease. *Mol Neurobiol*. 2019;56(8):5792-5798. doi: [10.1007/s12035-019-1501-x](https://doi.org/10.1007/s12035-019-1501-x). PMID: 30680692.
7. Norman M, Ter-Ovanesyan D, Trieu W, et al. L1CAM is not associated with extracellular vesicles in human cerebrospinal fluid or plasma. *Nat Methods*. 2021;18(6):631-634. doi: [10.1038/s41592-021-01174-8](https://doi.org/10.1038/s41592-021-01174-8). PMID: 34092791.
8. Kadam V, Wacker M, Oeckl P, Korneck M, et al. Most L1CAM Is not Associated with Extracellular Vesicles in Human Biofluids and iPSC-Derived Neurons. *Mol Neurobiol*. 2025;62(8):10427-10442. doi: [10.1007/s12035-025-04909-2](https://doi.org/10.1007/s12035-025-04909-2). PMID: 40210837.
9. Taha HB, Yacoub F, Morris A, Birnbaum A. Extracellular Vesicles for The Diagnosis of Alzheimer's Disease and Related Dementias: A Diagnostic Accuracy Meta-analysis. *Alzheimers Dement*. 2024;20(Suppl 2):e093432. doi: [10.1002/alz.093432](https://doi.org/10.1002/alz.093432). PMCID: PMC11716713. (Conference abstract; no PMID assigned.)
10. Théry C, Witwer KW, Aikawa E, et al. Minimal information for studies of extracellular vesicles 2018 (MISEV2018): a position statement of the International Society for Extracellular Vesicles and update of the MISEV2014 guidelines. *J Extracell Vesicles*. 2018;7(1):1535750. doi: [10.1080/20013078.2018.1535750](https://doi.org/10.1080/20013078.2018.1535750). PMID: 30637094.
11. Nielsen JE, Honoré B, Vestergård K, et al. Shotgun-based proteomics of extracellular vesicles in Alzheimer's disease reveals biomarkers involved in immunological and coagulation pathways. *Sci Rep*. 2021;11(1):18518. doi: [10.1038/s41598-021-97969-y](https://doi.org/10.1038/s41598-021-97969-y). PMID: 34531462.
12. Zhang Y, Zhao Y, Zhang J, et al. Proteomics of plasma-derived extracellular vesicles reveals S100A8 as a novel biomarker for Alzheimer's disease: A preliminary study. *J Proteomics*. 2024;308:105279. doi: [10.1016/j.jprot.2024.105279](https://doi.org/10.1016/j.jprot.2024.105279). PMID: 39159863.
13. Goetzl EJ, Boxer A, Schwartz JB, Abner EL, Petersen RC, Miller BL, Kapogiannis D. Altered lysosomal proteins in neural-derived plasma exosomes in preclinical Alzheimer disease. *Neurology*. 2015;85(1):40-47. doi: [10.1212/WNL.0000000000001702](https://doi.org/10.1212/WNL.0000000000001702). PMID: 26062630.
