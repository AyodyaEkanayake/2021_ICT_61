# Uncovering the Molecular Response of Rice to Drought Stress Using Differential Gene Expression Analysis

## 1. Project Overview

This research project investigates the molecular response of rice (*Oryza sativa*) to drought stress using differential gene expression analysis.

The main objective is to identify genes or microarray probes that show significant expression changes between control and drought-stressed rice samples.

## 2. Dataset Information

- **Dataset:** GSE115826
- **Source:** NCBI Gene Expression Omnibus (GEO)
- **Organism:** Oryza sativa (Rice)
- **Platform:** GPL6864
- **Data Type:** Microarray gene expression data

Dataset: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115826

## 3. Tools and Technologies

- R Programming Language
- RStudio
- GEOquery
- Biobase
- limma
- ggplot2
- pheatmap

## 4. Research Workflow

1. Download the dataset from NCBI GEO.
2. Extract expression values and sample metadata.
3. Identify control and drought-stressed samples.
4. Perform data quality checks and preprocessing.
5. Conduct exploratory data analysis.
6. Perform differential expression analysis using limma.
7. Identify upregulated and downregulated probes.
8. Generate volcano plots and expression heatmaps.
9. Obtain platform annotation information.
10. Export the analysis results.

## 5. Project Structure

```text
Rice_Drought_ICA3/
│
├── data/
├── results/
├── figures/
├── reports/
├── analysis.R
└── README.md
```

## 6. Main Analysis Methods

- Data preprocessing
- Exploratory Data Analysis (EDA)
- Multidimensional Scaling (MDS)
- Linear modelling using limma
- Empirical Bayes analysis
- Benjamini-Hochberg multiple testing correction

## 7. Significance Criteria

The following thresholds are used to identify significant probes:

- Adjusted P-value < 0.05
- Absolute log2 Fold Change ≥ 1

## 8. Generated Outputs

- All differential expression results
- Significant probes
- Upregulated probes
- Downregulated probes
- Annotated results
- Volcano plot
- Expression heatmap
- Sample metadata

## 9. Reproducibility

Run the `analysis.R` script in RStudio after installing the required R packages.

The dataset is downloaded from NCBI GEO using GEOquery.

## 10. Author

Undergraduate Research Project – BSc (Hons) in Information Technology

University of Vavuniya
