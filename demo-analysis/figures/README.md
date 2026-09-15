# Figures

This tree is **structured**. Figures live in one of three shapes, and the tree never goes deeper than three levels:

- **Stage-1 cohort figures** → `metadata/<family>/`
- **Stage-3 QC figures** → `qc/<family>/`
- **Stage-4+ results figures** → `analysis/<family>/<label>/`

The **family** is the name of the plot/analysis template that produced the figure. Processing state (e.g. `raw`, `normalized`, `log2`) goes in the **file stem**, never in a directory name. A figure attached to a finding has a stem that **starts with the finding id**.
