# The Findings Workflow — companion repository

Companion code and materials for the publication:

> **The Findings Workflow: An agentic, collaborative, provenance-tracked framework for scientific data analysis with Claude Code**
>
> Michael Riffle, Alex Zelter, Michael J. MacCoss

This repository lets readers of the paper see the Findings Workflow **in action, end to
end** — not just described, but as a complete, real run. It contains two things: the full
**transcript** of the collaborative Claude Code session that drove a sample analysis, and
the **actual project directory** that session produced, left exactly as the workflow wrote
it — every finding, figure, script, cached result, and piece of provenance in place.

## About the Findings Workflow

The Findings Workflow is an open-source, free [Claude Code](https://claude.com/claude-code)
plugin that addresses two recurring problems scientists face once the data are collected and
they are staring at a matrix of numbers: **confidence** (Is the data being normalized, batch
corrected, and tested correctly — especially where there is no canonical, boilerplate
solution?) and **provenance** (Exactly how was a result reached, and can it be reproduced?).
Using an LLM ad hoc can worsen both: convincing but flawed analyses, and a long chat history
that is not a reproducible method.

The workflow is opinionated and fully agentic. It delegates work such as statistics and
programming to multi-agent subsystems, **each with discrete, independent verification steps**,
and it **records findings automatically as they emerge**, complete with the provenance needed
to regenerate them (the script, the data version, the cached result). See the paper for the
full motivation and design.

The plugin lives at **[github.com/mriffle/findings-ai-collab-workflow](https://github.com/mriffle/findings-ai-collab-workflow)**.
Install it from within Claude Code:

```text
/plugin marketplace add mriffle/findings-ai-collab-workflow
/plugin install findings-workflow@findings-workflow
```

## What's in this repository

```
manuscript-findings-workflow/
├── README.md                            ← this file
├── findings-workflow-transcript-part1.txt   ┐ the full collaborative session transcript
├── findings-workflow-transcript-part2.txt   ┘ (human ⇄ agent), in two parts
├── figures/                             ← manuscript figure source files (Affinity Designer)
│   ├── figure1.af
│   └── figure2.af
└── demo-analysis/                       ← THE ACTUAL FINDINGS WORKFLOW RUN (see guide below)
```

- **`findings-workflow-transcript-part1.txt` / `-part2.txt`** — the complete record of the
  Claude Code session that carried out the sample analysis, from the first `init` through the
  exploration loop. This is the *collaborative* half of the paper's thesis: you can read the
  entire scientist ⇄ agent dialogue, the questions the workflow asked at each decision point,
  the subagents it dispatched, and the verification steps they ran.
- **`figures/`** — Affinity Designer (`.af`) source files for the manuscript's figures.
- **`demo-analysis/`** — the working directory in which the workflow ran. This is the
  *provenance* half of the thesis: everything the run produced, browsable in place. Guide
  below.

## The sample analysis

The demonstration dataset is a **pilot plasma extracellular-vesicle (EV) proteomics** study
(Mag-Net EV enrichment; data-independent-acquisition mass spectrometry): **40 individuals,
10 each** across four clinical groups — **ADD** (Alzheimer's disease dementia), **PDD**
(Parkinson's disease dementia), **PDCN** (Parkinson's, cognitively normal), and **HCN**
(healthy, cognitively normal). The scientific goal was to look for protein markers specific
to **ADD** versus other neurodegenerative groups.

This is a **demonstration of the workflow, not a scientific claim about the disease.** In
fact, one of the clearest things the run illustrates is the workflow's skepticism machinery
working as intended: the very first findings (0001–0005) record that diagnosis is
**confounded with collection site** and the data are small-*n* and already batch-corrected,
and every downstream discovery is therefore marked *exploratory / candidate* and
**confound-bounded** — never overstated. Reading the findings shows how the workflow keeps a
result honest about what it can and cannot prove.

## Guide to `demo-analysis/` — where everything is

The directory follows the Findings Workflow's standard layout. Each subdirectory has a
defined purpose, and the whole thing is designed to be read as a self-documenting record.

| Path | What it holds |
|------|---------------|
| **`CLAUDE.md`** | The workflow's standing instructions for this project — the rules the agent operated under the whole run. A good place to understand *how* the workflow behaves. |
| **`state/`** | Canonical, version-stamped project state: **`PROJECT.md`** (the science / study design), **`METADATA.md`** (verified column meanings, confounds), **`DATA_DESCRIPTION.md`** (matrix orientation, normalization/batch state, missingness), `color_registry.json`, and **`workflow.json`** (pipeline position + integrity-gate record). |
| **`findings/`** | The **findings graph** — the heart of the run. Thirteen findings (**`0001`–`0013`**), each a self-contained, provenance-tracked, evidence-bearing result with figures embedded inline. **`manifest.md`** is the index; **`exploration-log.md`** records what was tried. Start here: `findings/manifest.md`. |
| **`results/`** | Cached, **regenerable** analysis results — one directory per analysis (differential abundance, classification, Boruta, enrichment, GSEA, four-arm profiling, …), keyed by an input fingerprint, with a `manifest.md` registry. Findings point at these; nothing here is hand-edited. |
| **`figures/`** | Generated, publication-ready figures (vector **SVG** + 300 DPI **PNG** + a separate legend image), organized `metadata/`, `qc/`, and `analysis/`. Every figure is regenerable from a script, never hand-drawn. |
| **`scripts/`** | The analysis code: **`scratch/`** (exploratory scripts), **`promoted/`** (reviewed, tested scripts a finding may cite), **`figures/`** (figure-producing scripts), **`common/`** (shared modules), and `manifest.md`. |
| **`research/`** | Grounded **research-findings** — literature/biology background with independently **verified citations** (DOI/PMID), fact-checked before entering the corpus — plus `manifest.md`. |
| **`reports/`** | Generated reports (e.g. `stage3-qc-report.md`, the data-integrity QC report). |
| **`data/`** | The **raw input data** (read-only): the protein abundance matrix, sample metadata, and a `README.md`. |
| **`pyproject.toml`, `uv.lock`, `.python-version`** | The pinned, project-local Python environment used for every analysis (part of the provenance). `.venv/` and tool caches are git-ignored and regenerable. |

## Suggested reading paths

- **New to the workflow?** Read the paper, then `demo-analysis/CLAUDE.md` to see the rules
  the agent followed.
- **Want the science story?** `demo-analysis/state/PROJECT.md` → `demo-analysis/findings/manifest.md`
  → open individual findings (each has a plain-language summary, verdict, evidence, embedded
  figures, and caveats).
- **Want to see the human ⇄ agent collaboration?** Read the transcript files
  (`findings-workflow-transcript-part1.txt`, then `-part2.txt`).
- **Want to check reproducibility / provenance?** Pick any finding, follow its provenance
  block to the producing script (`demo-analysis/scripts/`), its cached result
  (`demo-analysis/results/`), and the pinned environment.

## How the pieces fit together

The paper's two claims map directly onto this repo. **Collaborative:** the transcript is the
complete dialogue that produced the work. **Provenance-tracked:** `demo-analysis/` is the
durable artifact — a finding is not a chat message but a script + a data version + a cached
result, regenerable by anyone later, with figures and caveats attached. You can start from
any figure or claim in a finding and trace it all the way back to the raw data.
