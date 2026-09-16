---
# Research-finding document — external knowledge, structured like a finding but for the literature/tools.
# Spec: doc 04.4. Stored in research/<slug>.md. A research finding without verified references is NOT accepted.

topic: "Is complement and/or coagulation-cascade activation a reported, replicated feature of Alzheimer's disease (AD) in blood/plasma and/or CSF? Evidence for mechanistic plausibility (complement-mediated synaptic pruning, AD mouse models), AD genetics (CLU, CR1 GWAS loci), complement proteins as blood/CSF biomarkers (direction and replication status), and coagulation-cascade involvement (fibrinogen, factor XII contact activation, fibrin deposition)."
type: research-publications
slug: complement-coagulation-alzheimers-blood-csf-biology
created: 2026-09-15
updated: 2026-09-15

status: reviewed
reviewed_by: research-reviewer (claude@ogdb.com)
reviewed_date: 2026-09-15

# Normalized domain-entity references this research is about (canonical IDs; gene symbol as label).
entities:
  - { db: uniprot, id: "P02745", label: "C1QA" }
  - { db: uniprot, id: "P02746", label: "C1QB" }
  - { db: uniprot, id: "P02747", label: "C1QC" }
  - { db: uniprot, id: "P09871", label: "C1S" }
  - { db: uniprot, id: "P01024", label: "C3" }
  - { db: uniprot, id: "P0C0L4", label: "C4A" }
  - { db: uniprot, id: "P17927", label: "CR1" }
  - { db: uniprot, id: "P11215", label: "ITGAM" }
  - { db: uniprot, id: "P10909", label: "CLU" }
  - { db: uniprot, id: "P08603", label: "CFH" }
  - { db: uniprot, id: "P02671", label: "FGA" }
  - { db: uniprot, id: "P02675", label: "FGB" }
  - { db: uniprot, id: "P02679", label: "FGG" }
  - { db: uniprot, id: "P00748", label: "F12" }
  - { db: uniprot, id: "P01042", label: "KNG1" }

# (required, non-empty) Every external claim traces to a reference here. Each is fact-checked by the
# research-reviewer: it must EXIST and SUPPORT the claim attributed to it. Hallucinated citations are
# the known failure mode this section + review exist to stop.
references:
  - id: "pmid:18083105"
    type: doi
    claim: "C1q, the initiating protein of the classical complement cascade, is expressed by postnatal neurons and localizes to synapses; mice deficient in C1q or the downstream complement protein C3 show large, sustained defects in developmental CNS synapse elimination, establishing complement as a molecular tag for synapse removal and proposing that this mechanism 'may become aberrantly reactivated in neurodegenerative disease.'"
    verified: true
    verified_by: research-reviewer (claude@ogdb.com)
  - id: "pmid:27033548"
    type: doi
    claim: "In AD mouse models, C1q is increased and associates with synapses before overt amyloid plaque deposition, and inhibition of C1q, C3, or the microglial complement receptor CR3 reduces microglial engulfment of synaptic material and early synapse loss, demonstrating that complement-microglia-mediated synaptic pruning is reactivated early in AD pathogenesis, not only as a late inflammatory epiphenomenon."
    verified: true
    verified_by: research-reviewer (claude@ogdb.com)
  - id: "pmid:28566429"
    type: doi
    claim: "In aged, plaque-rich APP/PS1 transgenic mice, genetic deletion of complement C3 protected against age-dependent hippocampal synapse and neuron loss and improved learning/memory performance despite the C3-deficient mice having MORE cerebral amyloid plaques, indicating that C3-dependent complement activation drives synaptic and neuronal loss downstream of, and independent of, plaque burden itself."
    verified: true
    verified_by: research-reviewer (claude@ogdb.com)
  - id: "pmid:18562603"
    type: doi
    claim: "In a different amyloid-precursor-protein (APP) transgenic AD mouse model, C3 deficiency instead led to up to twofold increased amyloid plaque burden, increased insoluble Abeta42, and significant hippocampal neuron loss at 17 months of age, indicating a beneficial, plaque-clearing role for C3 in this model — a result opposite in direction to the protective effect of C3 deletion reported in aged APP/PS1 mice, showing that complement's net effect on AD-related neuropathology is model-, age-, and context-dependent rather than uniformly deleterious."
    verified: true
    verified_by: research-reviewer (claude@ogdb.com)
  - id: "pmid:19734902"
    type: doi
    claim: "A two-stage genome-wide association study (GWAS) of over 16,000 individuals found genome-wide-significant association between Alzheimer's disease and a variant at the CLU (clusterin/APOJ) gene (rs11136000, combined P = 8.5 x 10^-10, odds ratio 0.86), alongside a locus near PICALM, establishing CLU as a novel AD risk gene."
    verified: true
    verified_by: research-reviewer (claude@ogdb.com)
  - id: "pmid:19734903"
    type: doi
    claim: "An independent, contemporaneous genome-wide association study of AD (2,032 French cases/5,328 controls, replicated in Belgian/Finnish/Italian/Spanish samples totaling 3,978 cases/3,297 controls) identified genome-wide-significant association with variants at CLU (rs11136000, OR = 0.86, P = 7.5 x 10^-9) and at CR1 (rs6656401, OR = 1.21, P = 3.7 x 10^-9) — CR1 encoding the complement component (3b/4b) receptor 1 — proposing roles for both genes in amyloid-beta clearance."
    verified: true
    verified_by: research-reviewer (claude@ogdb.com)
  - id: "pmid:20460622"
    type: doi
    claim: "A three-stage GWAS of more than 35,000 individuals states in its own framing that 'genome-wide association studies (GWAS) have recently identified CLU, PICALM, and CR1 as novel genes for late-onset Alzheimer disease,' and independently confirms the CLU and PICALM associations (P < .05) in an independent Spanish replication sample, providing early cross-cohort corroboration of the CLU/CR1/PICALM AD risk loci beyond their original discovery cohorts."
    verified: true
    verified_by: research-reviewer (claude@ogdb.com)
  - id: "pmid:24162737"
    type: doi
    claim: "A two-stage GWAS meta-analysis of 74,046 individuals of European ancestry (the International Genomics of Alzheimer's Project, IGAP) found 19 genome-wide-significant AD risk loci in total (including APOE), of which 11 were newly identified, substantially expanding the genome-wide-significant AD risk-locus roster beyond the original CLU/CR1/PICALM discovery in a much larger, independent sample."
    verified: true
    verified_by: research-reviewer (claude@ogdb.com)
  - id: "pmid:30820047"
    type: doi
    claim: "A genome-wide association meta-analysis of clinically diagnosed late-onset Alzheimer's disease (LOAD) in 94,437 individuals confirmed 20 previously reported LOAD risk loci, identified 5 new genome-wide-significant loci, and found that pathway analysis of AD risk genes implicates immunity, lipid metabolism, tau-binding proteins, and APP/amyloid-beta processing — i.e., immune-pathway genes (a category that includes complement genes such as CR1 and CLU established in earlier GWAS) remain enriched among AD risk loci in the largest case-control AD genetic studies."
    verified: true
    verified_by: research-reviewer (claude@ogdb.com)
  - id: "pmid:24794147"
    type: doi
    claim: "CR1 (complement component receptor 1), encoding a type-I transmembrane glycoprotein, is reviewed as 'one of the most important risk genes for late-onset Alzheimer's disease (LOAD),' with CR1 polymorphisms related to LOAD risk and CR1 proposed to affect AD pathogenesis via amyloid-beta pathology, tauopathy, immune dysfunction, and glial-mediated neuroinflammation."
    verified: true
    verified_by: research-reviewer (claude@ogdb.com)
  - id: "pmid:37480051"
    type: doi
    claim: "In a large cohort (early-onset AD n=912, late-onset AD n=492, controls n=504), plasma clusterin and C1q were significantly INCREASED (p < 0.001) and soluble CR1 and complement factor H were significantly REDUCED (p < 0.01) in AD versus controls; C1q was the single most predictive complement biomarker (AUC 0.655 for late-onset AD); and AD-associated GWAS SNPs in CR1 (rs6656401), C1S (rs3919533), and CFH (rs6664877) significantly influenced plasma levels of their respective encoded proteins, directly linking AD complement genetic risk loci to measurable plasma complement protein changes."
    verified: true
    verified_by: research-reviewer (claude@ogdb.com)
  - id: "pmid:22488444"
    type: doi
    claim: "In cerebrospinal fluid (CSF) from 43 AD, 42 MCI-progressing-to-AD (MCI-AD), 42 stable-MCI, and 44 control subjects, CSF levels of complement components C3 and C4 were elevated in AD compared with non-progressing MCI, and CSF CR1 was elevated in the combined MCI-AD/AD group, supporting 'aberrant complement regulation as a part in the AD process,' though the authors note these complement changes were not by themselves diagnostically useful."
    verified: true
    verified_by: research-reviewer (claude@ogdb.com)
  - id: "pmid:31628417"
    type: doi
    claim: "A systematic review and random-effects meta-analysis of 86 studies measuring CSF and/or peripheral blood complement/complement-regulator protein concentrations in AD versus healthy controls found that, in CSF, clusterin (SMD = 0.53, p < 0.005, I2 < 0.005%) and complement C3 (SMD = 0.45, p < 0.005, I2 = 68.40%) were significantly and fairly consistently higher in AD, while CSF C1q, CRP, SAP, and factor H were not significantly different; in peripheral blood, only CRP was significantly elevated in AD (SMD = 0.44, I2 = 93.81% — very high heterogeneity), while blood C3, C4, C1-inhibitor, SAP, factor H, and clusterin were NOT significantly different between AD and controls and were inconsistent between studies; the authors conclude that 'complement activity related to an AD diagnosis was not reflected consistently by the peripheral blood proteins investigated,' with the CSF clusterin and (less consistently) CSF C3 signals being the best-characterized replicated complement changes in AD."
    verified: true
    verified_by: research-reviewer (claude@ogdb.com)
  - id: "pmid:27567854"
    type: doi
    claim: "Using a multiplex assay to measure five complement proteins and four complement activation products in plasma from MCI, AD, and control donors, only one complement analyte (clusterin) differed significantly between AD and control plasma (controls 295 mg/L vs. AD 388 mg/L, p < 10^-5); three analytes (clusterin, factor I, terminal complement complex) differed significantly between MCI individuals who converted to dementia within a year versus non-converters; the paper explicitly notes that although 'numerous studies have suggested that measurement of individual complement proteins or activation products in cerebrospinal fluid or plasma is useful in diagnosis, prediction, or stratification, ... few have been replicated.'"
    verified: true
    verified_by: research-reviewer (claude@ogdb.com)
  - id: "pmid:25589717"
    type: doi
    claim: "In a case-control study comparing plasma complement factor H (CFH) levels in 317 AD cases and 254 controls (immunodiffusion assay, 80% power to detect a 23 mg/L between-group difference), no significant difference in plasma CFH was found between AD cases and controls, leading the authors to conclude that 'plasma CFH may not be a suitable biomarker for AD' — a directly negative/non-replicating result for an individual complement regulator in blood."
    verified: true
    verified_by: research-reviewer (claude@ogdb.com)
  - id: "pmid:20603455"
    type: doi
    claim: "In a combined proteomic/neuroimaging discovery-and-replication study (AddNeuroMed and the Baltimore Longitudinal Study of Aging), increased plasma clusterin/apolipoprotein J concentration was associated with entorhinal cortex atrophy, greater baseline AD severity, more rapid clinical progression, and predicted greater fibrillar amyloid-beta burden in the medial temporal lobe; APP/PS1 transgenic mice also showed increased plasma clusterin with age-dependent increases in brain clusterin and amyloid/clusterin colocalization in plaques."
    verified: true
    verified_by: research-reviewer (claude@ogdb.com)
  - id: "pmid:21467285"
    type: doi
    claim: "In a case-cohort study nested within the population-based Rotterdam Study, higher plasma clusterin was associated with increased likelihood of PREVALENT AD (odds ratio 1.63 per SD increase, 95% CI 1.21-2.20) and with greater disease severity (lower MMSE), but plasma clusterin levels were NOT associated with the risk of INCIDENT AD during prospective follow-up (adjusted HR 1.00, 95% CI 0.85-1.17), indicating that plasma clusterin tracks AD presence/severity cross-sectionally better than it predicts future AD risk."
    verified: true
    verified_by: research-reviewer (claude@ogdb.com)
  - id: "pmid:21824521"
    type: doi
    claim: "In a longitudinal MRI study of non-demented and MCI individuals, higher baseline plasma clusterin concentration was associated with slower rates of brain atrophy in the MCI group (but not in cognitively normal controls) across several brain regions, and in an autopsy sample, antemortem plasma clusterin concentration correlated significantly with postmortem clusterin concentration in the superior temporal gyrus, linking peripheral clusterin levels to brain-region-specific AD-relevant pathology."
    verified: true
    verified_by: research-reviewer (claude@ogdb.com)
  - id: "pmid:39455500"
    type: doi
    claim: "A 2025 review states that 'roles of complement in Alzheimer's disease (AD) and other dementias have emerged in recent years, supported by genetic, biomarker and pathological evidence and animal model studies,' summarizing complement's roles in brain homeostasis and disease and evaluating complement as a candidate AD therapeutic target."
    verified: true
    verified_by: research-reviewer (claude@ogdb.com)
  - id: "pmid:25775543"
    type: doi
    claim: "Plasma from human AD patients showed increased activation of coagulation factor XII (FXII), increased kallikrein activity, and increased cleavage of high-molecular-weight kininogen (HK) compared with controls; increased contact-system activation was also observed in AD mouse model plasma and in plasma from wild-type mice injected intravenously with amyloid-beta 42, demonstrating that circulating amyloid-beta can trigger the plasma FXII-driven contact-activation cascade (upstream of the intrinsic coagulation pathway) in human AD plasma."
    verified: true
    verified_by: research-reviewer (claude@ogdb.com)
  - id: "pmid:20547128"
    type: doi
    claim: "Fibrin clots formed in the presence of amyloid-beta are structurally abnormal and resistant to degradation; fibrin(ogen) was observed colocalizing with amyloid-positive blood vessels in both mouse and human AD brain samples, intravital imaging revealed abnormal thrombosis and fibrinolysis in AD mice, and depleting fibrinogen lessened cerebral amyloid angiopathy pathology and reduced cognitive impairment in AD mice, implicating fibrin(ogen) as a mechanistic contributor to AD."
    verified: true
    verified_by: research-reviewer (claude@ogdb.com)
  - id: "pmid:22869464"
    type: doi
    claim: "A review of amyloid-beta/fibrinogen interactions in AD states that 'a prothrombotic state in AD is evidenced by increased clot formation, decreased fibrinolysis, and elevated levels of coagulation factors and activated platelets,' and that fibrinogen abnormally deposits with amyloid-beta in brain parenchyma and cerebral blood vessels in AD, potentially contributing to amyloid deposition, reduced cerebral blood flow, and neuroinflammation."
    verified: true
    verified_by: research-reviewer (claude@ogdb.com)
  - id: "pmid:16269641"
    type: doi
    claim: "In the prospective, population-based Rotterdam Study (fibrinogen measured in 2,835 persons, 395 incident dementia cases over a mean 5.7-year follow-up), higher plasma fibrinogen was associated with increased risk of incident Alzheimer's disease (adjusted hazard ratio 1.25 per SD increase, 95% CI 1.04-1.49) and an even stronger association with vascular dementia (HR 1.76, 95% CI 1.34-2.30), whereas C-reactive protein was NOT associated with dementia risk, suggesting the fibrinogen-AD association reflects fibrinogen's hemostatic rather than purely inflammatory properties."
    verified: true
    verified_by: research-reviewer (claude@ogdb.com)
  - id: "pmid:18657226"
    type: doi
    claim: "Human AD brain tissue (but not non-demented brain tissue) showed diffuse fibrinogen and immunoglobulin (IgG) immunoreactivity indicative of blood-brain-barrier leakiness, with fibrinogen immunoreactivity appearing in association with amyloid-beta deposits and areas of microglial reactivity; in a rat model, intra-hippocampal amyloid-beta(1-42) injection caused time-dependent increases in fibrinogen co-staining with microgliosis, and defibrinogenation (ancrod) or microglial blockade (anti-Mac-1/CD11b) reduced amyloid-beta-plus-fibrinogen-driven inflammatory reactivity and was neuroprotective, mechanistically linking extravasated coagulation protein to AD-associated neuroinflammation and neuronal damage."
    verified: true
    verified_by: research-reviewer (claude@ogdb.com)
  - id: "pmid:34531462"
    type: doi
    claim: "In a separate small (n=10 AD/10 MCI/9 control) untargeted plasma-EV proteomics study (surveyed in depth in the companion research document ad-plasma-ev-proteomic-markers.md), coagulation factor XIII subunits F13A1 and F13B were among the most strongly upregulated proteins in AD-derived plasma EVs, alongside acute-phase/inflammatory and complement proteins — a single-study, EV-specific observation noted here only as a point of contact with (not independent replication of) the broader plasma/CSF complement-and-coagulation literature surveyed in this document."
    verified: true
    verified_by: research-reviewer (claude@ogdb.com)
---

# Research: Complement and coagulation-cascade activation as reported, replicated features of Alzheimer's disease biology in blood/plasma and CSF

## Summary

Complement activation is a genuine, mechanistically grounded, and genetically supported feature of Alzheimer's disease (AD) — not merely a laboratory curiosity of blood collection. Mouse-model work established that the classical complement proteins C1q and C3, acting through the microglial receptor CR3, mediate developmental synaptic pruning [pmid:18083105], that this same complement-microglia pathway is inappropriately reactivated early in AD pathogenesis to drive synapse loss before plaques accumulate [pmid:27033548], and that C3 signaling contributes to synaptic/neuronal loss in aged, plaque-rich mice even though the net effect of removing C3 is model- and age-dependent (protective in one APP/PS1 model [pmid:28566429], detrimental in an earlier APP-only model [pmid:18562603]). Independently, two of the earliest and most robustly replicated AD GWAS loci outside APOE are complement genes — CLU (clusterin) and CR1 (complement receptor 1) — identified in 2009 [pmid:19734902][pmid:19734903] and independently corroborated shortly after [pmid:20460622], with immune-pathway enrichment persisting in every subsequent, much larger AD GWAS meta-analysis [pmid:24162737][pmid:30820047], and a 2023 study directly linking AD-associated CR1/C1S/CFH genetic variants to altered levels of their encoded proteins in AD plasma [pmid:37480051]. At the biomarker level, the picture is real but heterogeneous, and this document is explicit about that heterogeneity rather than overstating it: a 2021 systematic review/meta-analysis of 86 studies found the complement signal is clearest and most consistently replicated in CSF (clusterin and, less consistently, C3 elevated in AD) but "not reflected consistently" across peripheral-blood complement proteins as a group, with CRP being the one blood protein that replicated well [pmid:31628417]. Individual-study results in blood are correspondingly mixed — plasma clusterin, C1q, soluble CR1, and factor H were all significantly altered (in disease-consistent directions) in one very large (n≈1,900) recent cohort [pmid:37480051] and clusterin alone differed in a smaller multiplex study [pmid:27567854], whereas plasma complement factor H showed no difference in a dedicated case-control study [pmid:25589717], and plasma clusterin tracked AD severity/state but not future incident risk in a prospective population cohort [pmid:21467285]. Coagulation-cascade involvement in AD is likewise supported by converging, largely independent lines of evidence: elevated plasma fibrinogen prospectively predicts incident AD risk in a population cohort [pmid:16269641], amyloid-beta directly triggers plasma factor XII-driven contact-system activation in human AD patient plasma [pmid:25775543], amyloid-beta/fibrin(ogen) interactions produce abnormal, degradation-resistant clots that colocalize with cerebral amyloid angiopathy [pmid:20547128][pmid:22869464], and fibrinogen extravasation into AD brain parenchyma colocalizes with amyloid deposits and reactive microglia [pmid:18657226]. Taken together, the literature supports treating an AD-associated, blood-detectable complement-and-coagulation signal as biologically plausible and consistent with an established (if still-heterogeneous, especially in blood) disease process, rather than something that should be dismissed a priori as a pre-analytical artifact.

## Detailed findings

### 1. Complement-mediated synaptic pruning: mechanistic plausibility from AD mouse models

The mechanistic case for complement in AD begins with normal brain development. Stevens et al. (2007, *Cell*) showed that C1q, the initiating protein of the classical complement cascade, is expressed by postnatal neurons and localizes to synapses, and that mice deficient in C1q or the downstream complement protein C3 fail to properly eliminate excess synaptic connections during development — establishing complement as a molecular "tag" for synapse elimination, and explicitly proposing that "complement-mediated synapse elimination may become aberrantly reactivated in neurodegenerative disease" [pmid:18083105]. This prediction was borne out for AD specifically by Hong et al. (2016, *Science*), who showed in AD mouse models that C1q increases and associates with synapses *before* overt amyloid plaque deposition, and that inhibiting C1q, C3, or the microglial complement receptor CR3 reduces microglial phagocytosis of synaptic material and blocks early synapse loss — establishing that complement-microglia synaptic pruning is reactivated at the earliest stages of AD pathogenesis, not just as a late-stage inflammatory bystander effect [pmid:27033548].

At later, plaque-rich disease stages, Shi et al. (2017, *Science Translational Medicine*) found that genetically deleting C3 in aged APP/PS1 mice protected against age-dependent synapse and neuron loss and improved cognitive performance, despite these C3-deficient mice having *more* cerebral amyloid plaques than their C3-intact counterparts — indicating that C3-dependent complement activation drives synaptic/neuronal damage downstream of, and at least partly independent of, raw plaque burden [pmid:28566429].

**This mechanism is not universally protective-to-block, and that complexity should be stated plainly.** An earlier study by Maier et al. (2008, *J Neurosci*) found the opposite net effect in a different APP transgenic model: C3 deficiency there led to *increased* (up to twofold) plaque burden and hippocampal neuron loss at 17 months, consistent with a beneficial, plaque-clearing role for C3 in that context [pmid:18562603]. The two mouse studies used different transgenic models (APP/PS1 vs. APP-only) and disease stages, and the discrepancy is a genuine, acknowledged feature of the field — complement's net effect on AD neuropathology is model-, age-, and context-dependent, not simply "more complement activation = worse pathology" in every setting. What is *not* in dispute across these studies is that complement (particularly C1q, C3, and CR3-bearing microglia) is mechanistically engaged with synapses and plaques in AD models; the open question the field continues to work on is the net balance of complement's clearance-promoting versus synapse/neuron-damaging effects at different disease stages.

### 2. AD genetics: CLU and CR1 as replicated complement-pathway GWAS loci

Independent of the mouse mechanistic work, human genetics converged on complement genes as AD risk factors starting in 2009. Two back-to-back genome-wide association studies (GWAS) each independently identified CLU (clusterin/APOJ) as a genome-wide-significant AD risk locus [pmid:19734902][pmid:19734903], and Lambert et al. additionally identified CR1 — which encodes the complement component (3b/4b) receptor 1 — as genome-wide significant (rs6656401, OR = 1.21, P = 3.7 x 10^-9), proposing that both CLU and CR1 act via roles in amyloid-beta clearance [pmid:19734903]. A third, larger, multi-stage GWAS (Seshadri et al. 2010, *JAMA*) explicitly frames CLU, PICALM, and CR1 as the set of "novel genes for late-onset Alzheimer disease" recently identified by GWAS, and independently replicates the CLU/PICALM associations in a further Spanish sample — early cross-cohort corroboration beyond the original discovery cohorts [pmid:20460622]. Every subsequent, substantially larger AD GWAS meta-analysis has continued to find immune-pathway genes enriched among genome-wide-significant AD risk loci: the 2013 IGAP meta-analysis of 74,046 individuals found 19 genome-wide-significant loci in total (11 newly identified) [pmid:24162737], and the 2019 Kunkle et al. meta-analysis of 94,437 individuals confirmed 20 previously reported LOAD risk loci and found pathway analysis implicating immunity alongside lipid metabolism and APP/amyloid-beta processing [pmid:30820047]. CR1 itself has been the subject of a dedicated review summarizing its role as "one of the most important risk genes" for late-onset AD, with proposed downstream effects on amyloid-beta pathology, tauopathy, immune dysfunction, and glial-mediated neuroinflammation [pmid:24794147].

Crucially for the question of whether AD complement genetics should be expected to produce a measurable *blood* signal, Veteleanu et al. (2023, *J Neuroinflammation*) directly closed this loop: in a large cohort (early-onset AD n=912, late-onset AD n=492, controls n=504), AD-associated GWAS SNPs in CR1 (rs6656401), C1S (rs3919533), and CFH (rs6664877) significantly influenced plasma levels of their respective encoded proteins, and the same cohort showed plasma clusterin and C1q significantly *increased* and soluble CR1 and factor H significantly *reduced* in AD versus controls (p < 0.001 and p < 0.01, respectively), with C1q the single most predictive complement biomarker (AUC 0.655 for late-onset AD) [pmid:37480051]. This is the most direct available evidence that AD complement genetic risk translates into measurable, disease-associated plasma complement protein changes, and it establishes plausible mechanistic direction for at least four specific complement proteins in plasma (clusterin up, C1q up, sCR1 down, CFH down).

### 3. Complement proteins as AD biomarkers in blood/plasma and CSF: what replicates and what does not

This is the section where honesty about mixed and heterogeneous evidence matters most, and the single most load-bearing reference here is a dedicated systematic review and meta-analysis. Krance et al. (2021, *Molecular Psychiatry*) systematically reviewed 2,966 records and meta-analyzed 86 original studies measuring CSF and/or peripheral blood complement or complement-regulator protein concentrations in AD versus healthy controls [pmid:31628417]. Their findings, stated plainly:

- **In CSF**, clusterin (standardized mean difference, SMD = 0.53, p < 0.005, very low heterogeneity I2 < 0.005%) and complement C3 (SMD = 0.45, p < 0.005, but higher heterogeneity I2 = 68.4%) were significantly and fairly consistently **higher** in AD than controls. CSF C1q, C-reactive protein (CRP), serum amyloid protein (SAP), and factor H were **not** significantly different.
- **In peripheral blood**, only CRP was significantly elevated in AD (SMD = 0.44), but with very high between-study heterogeneity (I2 = 93.81%). Blood C3, C4, C1-inhibitor, SAP, factor H, and clusterin were **not** significantly different between AD and controls and were described as inconsistent between studies.
- The authors' own conclusion: complement pathway activity is elevated in AD, "best characterized by increased CSF clusterin concentrations and less consistently by CSF C3 concentrations," while "complement activity related to an AD diagnosis was not reflected consistently by the peripheral blood proteins investigated" [pmid:31628417].

This meta-analytic picture is corroborated and refined by individual studies. In CSF, Daborg et al. (2012) found C3 and C4 elevated in AD compared with non-progressing MCI, and CR1 elevated in the combined MCI-progressing-to-AD/AD group, describing this as support for "aberrant complement regulation" in the AD disease process — while explicitly noting the changes were not by themselves diagnostically useful [pmid:22488444]. In plasma, results for individual analytes diverge across studies:

- **Clusterin** is the most consistently reported plasma complement-pathway protein: it was the only analyte among five complement proteins and four activation products that differed significantly between AD and control plasma in a targeted multiplex study (controls 295 mg/L vs. AD 388 mg/L, p < 10^-5) [pmid:27567854]; it was significantly elevated and predictive in the large Veteleanu et al. cohort [pmid:37480051]; it tracked entorhinal cortex atrophy, baseline severity, rapid clinical progression, and fibrillar amyloid burden in a combined proteomic/neuroimaging study with mouse-model corroboration [pmid:20603455]; and it was associated with brain atrophy rate in MCI and correlated with postmortem brain clusterin concentration in a separate cohort [pmid:21824521]. However, a population-based prospective study (the Rotterdam Study) found that while plasma clusterin was associated with *prevalent* AD (OR 1.63 per SD) and disease severity, it was **not** associated with the risk of *incident* AD during follow-up (adjusted HR 1.00) [pmid:21467285] — an important nuance: plasma clusterin appears to track AD state/severity better than it predicts future onset, and should not be treated as an established AD risk-prediction biomarker even though it reliably differs cross-sectionally.
- **Complement factor H (CFH)** is a clear example of a plasma complement marker that has **not** replicated: a dedicated case-control study (317 AD cases, 254 controls, adequately powered) found no significant difference in plasma CFH and concluded "plasma CFH may not be a suitable biomarker for AD" [pmid:25589717] — directly contradicting the (opposite-direction, reduced-in-AD) CFH finding from the much larger Veteleanu et al. cohort [pmid:37480051]. The Hakobyan et al. multiplex study also explicitly cautions that of the many complement biomarkers proposed in the literature, "few have been replicated" [pmid:27567854].
- **C1q**, **soluble CR1**, and **C3/C4** show significant, disease-consistent-direction changes in specific large studies (C1q and clusterin up, sCR1 down [pmid:37480051]; C3/C4/CR1 up in CSF [pmid:22488444]) but were not confirmed as consistently different in blood across the broader meta-analyzed literature [pmid:31628417].

**Bottom line for this subsection, stated without overstatement:** the complement-in-AD biomarker signal is real and replicated at the level of CSF clusterin (and, less consistently, CSF C3) and at the level of AD genetics-to-plasma-protein linkage for CR1/C1S/CFH variants [pmid:37480051]. It is markedly less consistent when peripheral blood complement proteins are meta-analyzed as a group [pmid:31628417], and specific plasma analytes (CFH in particular) have produced directly conflicting results across studies [pmid:25589717][pmid:37480051]. An AD-associated complement signal detected in a plasma-derived sample is therefore biologically plausible and has real precedent in the literature, but individual-protein-level claims should be treated as provisional pending independent replication, consistent with how this literature itself describes its own reproducibility problem.

### 4. Coagulation-cascade involvement in AD

Independent of complement, the blood coagulation cascade has its own body of evidence linking it to AD, spanning population epidemiology, human patient plasma biochemistry, and mechanistic mouse/tissue work.

**Prospective, population-based evidence:** In the Rotterdam Study, a prospective population-based cohort, higher plasma fibrinogen (measured in 2,835 persons, with 395 incident dementia cases over a mean 5.7-year follow-up) was associated with increased risk of incident Alzheimer's disease (adjusted HR 1.25 per SD increase, 95% CI 1.04-1.49) and an even stronger association with vascular dementia (HR 1.76) — while C-reactive protein, a general inflammatory marker, was *not* associated with dementia risk, leading the authors to argue the fibrinogen-AD association reflects fibrinogen's hemostatic rather than purely inflammatory role [pmid:16269641]. This is a genuinely prospective, population-based, replication-grade design (not merely a cross-sectional case-control comparison), which strengthens confidence in a real fibrinogen-AD association.

**Human AD plasma biochemistry:** Zamolodchikov et al. (2015, *PNAS*) directly measured contact-system activation in human AD patient plasma and found increased activation of coagulation factor XII (FXII), increased kallikrein activity, and increased cleavage of high-molecular-weight kininogen (HK) compared with controls; the same pattern was reproduced in AD mouse model plasma and in wild-type mice injected intravenously with amyloid-beta 42, demonstrating that circulating amyloid-beta can directly trigger the plasma FXII-driven contact-activation cascade — a cascade that lies immediately upstream of, and cross-talks with, the intrinsic coagulation pathway — in human AD plasma [pmid:25775543].

**Mechanistic amyloid-beta/fibrin(ogen) interaction:** Cortes-Canteli et al. (2010, *Neuron*) showed that fibrin clots formed in the presence of amyloid-beta are structurally abnormal and resistant to degradation in vitro and in vivo; fibrin(ogen) colocalized with amyloid-positive blood vessels in both mouse and human AD brain samples; intravital imaging revealed abnormal thrombosis and fibrinolysis in AD mice; and depleting fibrinogen reduced cerebral amyloid angiopathy pathology and improved cognitive performance in AD mice [pmid:20547128]. A subsequent review by the same group summarizes this body of work, stating plainly that "a prothrombotic state in AD is evidenced by increased clot formation, decreased fibrinolysis, and elevated levels of coagulation factors and activated platelets," with fibrinogen abnormally co-depositing with amyloid-beta in brain parenchyma and cerebral vessels [pmid:22869464].

**Fibrinogen extravasation and neuroinflammation in AD brain tissue:** Ryu and McLarnon (2009, *J Cell Mol Med*) found that human AD brain tissue — but not non-demented brain tissue — showed diffuse fibrinogen and immunoglobulin (IgG) immunoreactivity indicative of blood-brain-barrier leakiness, with fibrinogen colocalizing with amyloid-beta deposits and areas of microglial reactivity. In a rat model, intra-hippocampal amyloid-beta(1-42) injection produced time-dependent fibrinogen/microgliosis co-staining, and either defibrinogenation (the snake-venom-derived agent ancrod) or microglial blockade (anti-Mac-1/CD11b) reduced amyloid-beta-plus-fibrinogen-driven inflammatory and vascular damage and was neuroprotective — mechanistically tying extravasated coagulation protein directly to AD-associated neuroinflammation and neuronal damage [pmid:18657226].

Taken together, these four independent lines of evidence (population-epidemiological, direct human-plasma biochemical, in vitro/in vivo mechanistic, and human-brain-tissue immunohistochemical) converge on the same conclusion: coagulation-cascade activation, and specifically fibrinogen/fibrin and the FXII-driven contact system, is a real, multiply-supported feature of AD pathophysiology extending into blood/plasma, not an isolated or single-study finding.

### 5. A note on adjacent EV-proteomics evidence (not the focus of this document)

A separate research document in this project's corpus (`research/ad-plasma-ev-proteomic-markers.md`) surveys plasma extracellular-vesicle (EV) proteomic markers of AD in depth, including a small (n=10 AD/10 MCI/9 control) untargeted MS study of plasma EVs in which coagulation factor XIII subunits F13A1 and F13B were among the most strongly upregulated proteins in AD, alongside acute-phase/inflammatory and complement proteins [pmid:34531462]. That finding is a single-study, EV-specific observation and is not treated here as independent replication of the broader complement/coagulation biology surveyed above — it is noted only as a point of contact between the EV-proteomics literature and the blood/CSF complement-and-coagulation biology this document establishes.

### 6. Synthesis: is a plasma-EV "complement and coagulation cascades" AD signal biologically plausible?

Based on the literature verified above, the answer is yes, with calibrated confidence that varies by claim:

- **High confidence, well-replicated:** complement (C1q/C3/CR3) mechanistically mediates synaptic pruning and is reactivated early in AD mouse models [pmid:18083105][pmid:27033548]; CLU and CR1 are genome-wide-significant, independently replicated AD genetic risk loci, with immune-pathway enrichment persisting in every larger subsequent GWAS [pmid:19734902][pmid:19734903][pmid:20460622][pmid:24162737][pmid:30820047]; CSF clusterin is consistently elevated in AD by meta-analysis [pmid:31628417]; and elevated plasma fibrinogen prospectively predicts incident AD risk in a population cohort, independent of general inflammation (CRP) [pmid:16269641].
- **Moderate confidence, plausible but heterogeneous:** peripheral blood complement protein levels as a category are not consistently different between AD and controls by meta-analysis [pmid:31628417], but specific large, genetically-anchored cohort studies do find significant, direction-consistent plasma complement changes (clusterin and C1q up; sCR1 and CFH down) tied to AD risk-allele genotype [pmid:37480051], and amyloid-beta directly activates the plasma FXII contact/coagulation-adjacent system in AD patient plasma [pmid:25775543].
- **Explicitly contested or non-replicating at the single-analyte level:** plasma CFH shows directly conflicting results across studies (no difference in one adequately powered case-control study [pmid:25589717] versus significantly reduced in a much larger cohort [pmid:37480051]); plasma clusterin tracks AD state/severity but not incident risk [pmid:21467285]; and complement's net mechanistic effect on plaque burden/neurodegeneration is mouse-model-dependent, with opposite-direction results reported for C3 deficiency in different APP transgenic lines [pmid:28566429][pmid:18562603].

None of this literature makes a specific claim about DIA-MS-detected complement/coagulation proteins in Mag-Net-enriched plasma EVs specifically — that link is not established by any source found in this search and should not be inferred from it. What this document does establish is that an AD-associated activation of complement and/or coagulation pathways detectable in blood-derived material is a real, mechanistically grounded, genetically supported, and (for several specific proteins/pathways) independently replicated feature of AD biology — sufficient grounds that such a signal in a plasma-EV proteomics dataset should not be dismissed a priori as a pre-analytical or collection artifact, while individual-protein directionality claims should still be checked against the mixed evidence summarized above rather than assumed uniform.

## Methods / sources consulted

- **Primary-source discovery and verification via PubMed E-utilities** (`esearch`/`esummary`/`efetch`, NCBI NLM, `eutils.ncbi.nlm.nih.gov`, accessed via `curl` in this session). Every PMID cited above was located via targeted `esearch` queries (author + specific finding/title fragment), its existence and identity confirmed via `esummary` (title, journal, pubdate, first author), and its full abstract retrieved via `efetch` (`rettype=abstract&retmode=text`) and read in full before any claim was attributed to it. No PMID or claim in this document was taken from memory or from a search-result snippet alone; every `claim` field above closely paraphrases or quotes language actually present in the corresponding PubMed abstract.
- **Search strategy:** iterative queries covering four sub-areas per the task scope: (1) complement-mediated synaptic pruning and AD mouse-model mechanism (Stevens, Hong, Shi, Maier); (2) AD GWAS loci implicating complement genes (Harold, Lambert x2, Seshadri, Kunkle, plus the CR1 review by Zhu and the genotype-to-plasma-protein study by Veteleanu); (3) complement proteins as blood/plasma/CSF AD biomarkers, prioritizing systematic reviews/meta-analyses over individual studies to establish replication status (Krance meta-analysis; Daborg CSF; Hakobyan, Williams, Thambisetty x2, Schrijvers plasma); (4) coagulation-cascade involvement in AD (Zamolodchikov, Cortes-Canteli x2, van Oijen, Ryu & McLarnon).
- **UniProt accession verification:** every accession in the `entities` list above was fetched directly from the UniProt REST API (`https://rest.uniprot.org/uniprotkb/<ACCESSION>.txt`) in this session, and the returned `ID` (entry name) and `GN` (gene name) lines were read directly from the response to confirm the accession corresponds to the intended canonical human Swiss-Prot entry (e.g., confirmed P17927 = CR1_HUMAN/CR1, P10909 = CLUS_HUMAN/CLU, P08603 = CFAH_HUMAN/CFH, P00748 = FA12_HUMAN/F12, P01042 = KNG1_HUMAN/KNG1).
- **Deliberately excluded / not pursued:** a full-text check of the very large multi-author GWAS papers (Lambert 2013, Kunkle 2019, and an additional 2022 GWAS by Bellenguez et al., PMID 35379992) was limited to their abstracts; the Bellenguez 2022 GWAS abstract was reviewed but **not cited** in this document because its abstract text does not itself name CR1 or complement specifically (it highlights microglia and TNF-pathway involvement), so a complement-specific claim could not be verified from it without full-text access — flagged here as a claim intentionally left out rather than approximated.
- **Not found / not independently verified in this search:** a study directly measuring the specific complement activation fragments C3b or iC3b as blood/CSF AD biomarkers (the task's example markers) was not identified with a specific, checkable citation; Hakobyan et al. [pmid:27567854] measured unspecified "complement activation products" in plasma but the abstract does not name C3b/iC3b specifically, so no claim naming those two fragments by name is made in this document.

## References

1. Stevens B, Allen NJ, Vazquez LE, et al. The classical complement cascade mediates CNS synapse elimination. *Cell*. 2007;131(6):1164-1178. doi: [10.1016/j.cell.2007.10.036](https://doi.org/10.1016/j.cell.2007.10.036). PMID: 18083105.
2. Hong S, Beja-Glasser VF, Nfonoyim BM, et al. Complement and microglia mediate early synapse loss in Alzheimer mouse models. *Science*. 2016;352(6286):712-716. doi: [10.1126/science.aad8373](https://doi.org/10.1126/science.aad8373). PMID: 27033548.
3. Shi Q, Chowdhury S, Ma R, et al. Complement C3 deficiency protects against neurodegeneration in aged plaque-rich APP/PS1 mice. *Sci Transl Med*. 2017;9(392):eaaf6295. doi: [10.1126/scitranslmed.aaf6295](https://doi.org/10.1126/scitranslmed.aaf6295). PMID: 28566429.
4. Maier M, Peng Y, Jiang L, Seabrook TJ, Carroll MC, Lemere CA. Complement C3 deficiency leads to accelerated amyloid beta plaque deposition and neurodegeneration and modulation of the microglia/macrophage phenotype in amyloid precursor protein transgenic mice. *J Neurosci*. 2008;28(25):6333-6341. doi: [10.1523/JNEUROSCI.0829-08.2008](https://doi.org/10.1523/JNEUROSCI.0829-08.2008). PMID: 18562603.
5. Harold D, Abraham R, Hollingworth P, et al. Genome-wide association study identifies variants at CLU and PICALM associated with Alzheimer's disease. *Nat Genet*. 2009;41(10):1088-1093. doi: [10.1038/ng.440](https://doi.org/10.1038/ng.440). PMID: 19734902.
6. Lambert JC, Heath S, Even G, et al. Genome-wide association study identifies variants at CLU and CR1 associated with Alzheimer's disease. *Nat Genet*. 2009;41(10):1094-1099. doi: [10.1038/ng.439](https://doi.org/10.1038/ng.439). PMID: 19734903.
7. Seshadri S, Fitzpatrick AL, Ikram MA, et al. Genome-wide analysis of genetic loci associated with Alzheimer disease. *JAMA*. 2010;303(18):1832-1840. doi: [10.1001/jama.2010.574](https://doi.org/10.1001/jama.2010.574). PMID: 20460622.
8. Lambert JC, Ibrahim-Verbaas CA, Harold D, et al. Meta-analysis of 74,046 individuals identifies 11 new susceptibility loci for Alzheimer's disease. *Nat Genet*. 2013;45(12):1452-1458. doi: [10.1038/ng.2802](https://doi.org/10.1038/ng.2802). PMID: 24162737.
9. Kunkle BW, Grenier-Boley B, Sims R, et al. Genetic meta-analysis of diagnosed Alzheimer's disease identifies new risk loci and implicates Aβ, tau, immunity and lipid processing. *Nat Genet*. 2019;51(3):414-430. doi: [10.1038/s41588-019-0358-2](https://doi.org/10.1038/s41588-019-0358-2). PMID: 30820047.
10. Zhu XC, Yu JT, Jiang T, Wang P, Cao L, Tan L. CR1 in Alzheimer's disease. *Mol Neurobiol*. 2015;51(2):753-765. doi: [10.1007/s12035-014-8723-8](https://doi.org/10.1007/s12035-014-8723-8). PMID: 24794147.
11. Veteleanu A, Stevenson-Hoare J, Keat S, et al. Alzheimer's disease-associated complement gene variants influence plasma complement protein levels. *J Neuroinflammation*. 2023;20(1):169. doi: [10.1186/s12974-023-02850-6](https://doi.org/10.1186/s12974-023-02850-6). PMID: 37480051.
12. Daborg J, Andreasson U, Pekna M, et al. Cerebrospinal fluid levels of complement proteins C3, C4 and CR1 in Alzheimer's disease. *J Neural Transm (Vienna)*. 2012;119(7):789-797. doi: [10.1007/s00702-012-0797-8](https://doi.org/10.1007/s00702-012-0797-8). PMID: 22488444.
13. Krance SH, Wu CY, Zou Y, et al. The complement cascade in Alzheimer's disease: a systematic review and meta-analysis. *Mol Psychiatry*. 2021;26(10):5532-5541. doi: [10.1038/s41380-019-0536-8](https://doi.org/10.1038/s41380-019-0536-8). PMID: 31628417.
14. Hakobyan S, Harding K, Aiyaz M, et al. Complement Biomarkers as Predictors of Disease Progression in Alzheimer's Disease. *J Alzheimers Dis*. 2016;54(2):707-716. doi: [10.3233/JAD-160420](https://doi.org/10.3233/JAD-160420). PMID: 27567854.
15. Williams MA, Haughton D, Stevenson M, Craig D, Passmore AP, Silvestri G. Plasma Complement factor H in Alzheimer's Disease. *J Alzheimers Dis*. 2015;45(2):369-372. doi: [10.3233/JAD-142742](https://doi.org/10.3233/JAD-142742). PMID: 25589717.
16. Thambisetty M, Simmons A, Velayudhan L, et al. Association of plasma clusterin concentration with severity, pathology, and progression in Alzheimer disease. *Arch Gen Psychiatry*. 2010;67(7):739-748. doi: [10.1001/archgenpsychiatry.2010.78](https://doi.org/10.1001/archgenpsychiatry.2010.78). PMID: 20603455.
17. Schrijvers EM, Koudstaal PJ, Hofman A, Breteler MM. Plasma clusterin and the risk of Alzheimer disease. *JAMA*. 2011;305(13):1322-1326. doi: [10.1001/jama.2011.381](https://doi.org/10.1001/jama.2011.381). PMID: 21467285.
18. Thambisetty M, An Y, Kinsey A, et al. Plasma clusterin concentration is associated with longitudinal brain atrophy in mild cognitive impairment. *Neuroimage*. 2012;59(1):212-217. doi: [10.1016/j.neuroimage.2011.07.056](https://doi.org/10.1016/j.neuroimage.2011.07.056). PMID: 21824521.
19. Daskoulidou N, Carpanini SM, Zelek WM, Morgan BP. Involvement of Complement in Alzheimer's Disease: From Genetics Through Pathology to Therapeutic Strategies. *Curr Top Behav Neurosci*. 2025;69:3-24. doi: [10.1007/7854_2024_524](https://doi.org/10.1007/7854_2024_524). PMID: 39455500.
20. Zamolodchikov D, Chen ZL, Conti BA, Renné T, Strickland S. Activation of the factor XII-driven contact system in Alzheimer's disease patient and mouse model plasma. *Proc Natl Acad Sci U S A*. 2015;112(13):4068-4073. doi: [10.1073/pnas.1423764112](https://doi.org/10.1073/pnas.1423764112). PMID: 25775543.
21. Cortes-Canteli M, Paul J, Norris EH, et al. Fibrinogen and beta-amyloid association alters thrombosis and fibrinolysis: a possible contributing factor to Alzheimer's disease. *Neuron*. 2010;66(5):695-709. doi: [10.1016/j.neuron.2010.05.014](https://doi.org/10.1016/j.neuron.2010.05.014). PMID: 20547128.
22. Cortes-Canteli M, Zamolodchikov D, Ahn HJ, Strickland S, Norris EH. Fibrinogen and altered hemostasis in Alzheimer's disease. *J Alzheimers Dis*. 2012;32(3):599-608. doi: [10.3233/JAD-2012-120820](https://doi.org/10.3233/JAD-2012-120820). PMID: 22869464.
23. van Oijen M, Witteman JC, Hofman A, Koudstaal PJ, Breteler MM. Fibrinogen is associated with an increased risk of Alzheimer disease and vascular dementia. *Stroke*. 2005;36(12):2637-2641. doi: [10.1161/01.STR.0000189721.31432.26](https://doi.org/10.1161/01.STR.0000189721.31432.26). PMID: 16269641.
24. Ryu JK, McLarnon JG. A leaky blood-brain barrier, fibrinogen infiltration and microglial reactivity in inflamed Alzheimer's disease brain. *J Cell Mol Med*. 2009;13(9A):2911-2925. doi: [10.1111/j.1582-4934.2008.00434.x](https://doi.org/10.1111/j.1582-4934.2008.00434.x). PMID: 18657226.
25. Nielsen JE, Honoré B, Vestergård K, et al. Shotgun-based proteomics of extracellular vesicles in Alzheimer's disease reveals biomarkers involved in immunological and coagulation pathways. *Sci Rep*. 2021;11(1):18518. doi: [10.1038/s41598-021-97969-y](https://doi.org/10.1038/s41598-021-97969-y). PMID: 34531462. (Noted in passing per §5 above; surveyed in full in the companion document `research/ad-plasma-ev-proteomic-markers.md`.)
