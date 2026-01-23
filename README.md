# ML Bioinformatics Workflow for (CKD, CVD, Hypertension) Gene Expression Analysis

Scripts and datasets for gene-expression exploration across CVD, CKD, and hypertension: preprocessing, ANOVA and mRMR feature selection, overlap plotting, and pathway/heatmap visuals.

## Data

- Folder: 0_datasets/
- Shape: samples as rows, genes as columns, final column = target/label.
- Add datasets by dropping more CSVs in the same layout.

## Run

- ANOVA (interactive threshold or top-N): `python 2_anova-feature-selection/2_feature_selection_anova.py`
- mRMR (edit `dataset_name` in script): `python 3_mrmr-feature-selection/3_mRMR.py`
- Overlap plot: `python 4_common_percentage.py`

## Scripts at a glance

- 1_get-genes/: quick exploration and filtering helpers.
- 2_anova-feature-selection/: compute ANOVA scores and save selected-feature CSVs plus reports.
- 3_mrmr-feature-selection/: rank top K genes with mRMR.
- R figures: 5_PathWays.R, 6_Heatmap.R, 7_reactom_pathway.R, 8_go_analysis.R, 9_box_plot.R.

## Requirements

- Python 3.9+, packages: pandas, numpy, matplotlib, mrmr.
- R scripts assume common tidyverse/bioconductor plotting and pathway libs.

## Authors

Francis Rudra D Cruze<sup>1</sup> (francisrudra@gmail.com) <br>
Md. Faruk Hosen<sup>2</sup> <br>
Mohammad Badrul Alam Miah<sup>3</sup> <br>
A. S. M. Sanwar Hosen<sup>4</sup> <br>

Affiliations

1. Department of Computer Science and Engineering, East West University, Aftabnagar, Dhaka 1212, Bangladesh.
2. Department of Computer Science and Engineering, Netrokona University, Netrokona 2400, Bangladesh.
3. Department of Information and Communication Technology, Mawlana Bhashani Science and Technology University, Tangail 1902, Bangladesh.
4. Department of Artificial Intelligence and Big Data, Woosong University, Daejeon 34606, South Korea.
