
# Uncovering the Molecular Response of Rice to Drought Stress Using Differential Gene Expression Analysis

## 1. Project Overview

This project investigates the molecular response of rice (*Oryza sativa*) to drought stress using differential gene expression analysis.

The main objective is to identify genes or microarray probes that exhibit significant expression differences between control and drought-stressed rice samples.

The analysis uses publicly available microarray data obtained from the National Center for Biotechnology Information (NCBI) Gene Expression Omnibus (GEO).

## 2. Dataset Information

| Attribute | Description |
|---|---|
| Dataset ID | GSE115826 |
| Source | NCBI Gene Expression Omnibus (GEO) |
| Organism | *Oryza sativa* (Rice) |
| Platform | GPL6864 |
| Data Type | Microarray gene expression data |
| Experimental Conditions | Control and Drought Stress |

**Dataset Link:**  
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115826

## 3. Tools and Technologies

The following tools and R packages are used in this project:

- **R** – Statistical computing and data analysis
- **RStudio** – Integrated development environment
- **GEOquery** – Retrieving gene expression data from NCBI GEO
- **Biobase** – Managing and processing gene expression datasets
- **limma** – Differential expression analysis
- **ggplot2** – Data visualization
- **pheatmap** – Generating expression heatmaps

## 4. Research Workflow

The project follows these main analytical stages:

1. Retrieve the GSE115826 dataset from NCBI GEO.
2. Extract gene expression values and sample metadata.
3. Identify control and drought-stressed samples.
4. Perform data quality checks and preprocessing.
5. Conduct exploratory data analysis (EDA).
6. Explore sample relationships using Multidimensional Scaling (MDS).
7. Perform differential expression analysis using the limma package.
8. Apply empirical Bayes moderation and multiple-testing correction.
9. Identify significant, upregulated, and downregulated probes.
10. Retrieve platform annotation information.
11. Generate visualizations, including volcano plots and expression heatmaps.
12. Export the analysis results and prepare the final report.

## 5. Project Structure

```text
Rice_Drought_ICA3/
│
├── code/
│   └── analysis.R
│
├── data/
│
├── figures/
│
├── results/
│
├── reports/
│   ├── final_report.Rmd
│   └── final_report.html
│
├── README.md
│
└── Rice_Drought_ICA3.Rproj
```

## 6. Main Analysis Methods

The following methods are applied throughout the analysis:

- Data preprocessing and quality assessment
- Exploratory Data Analysis (EDA)
- Descriptive statistical analysis
- Multidimensional Scaling (MDS)
- Linear modelling using the limma package
- Empirical Bayes moderation
- Benjamini-Hochberg multiple-testing correction
- Probe annotation and interpretation
- Data visualization

## 7. Differential Expression Significance Criteria

Probes are classified as significant based on the following criteria:

- **Adjusted P-value < 0.05**
- **Absolute log fold change (|logFC|) ≥ 1**

The adjusted P-value threshold is used to identify statistically significant expression differences while accounting for multiple hypothesis testing.

The original expression measurements deposited in the GEO Series Matrix file are retained during preprocessing. These measurements correspond to Agilent Feature Extraction LogRatio values. No additional logarithmic scale conversion is applied.

The exact logarithm base of the deposited LogRatio values is not assumed without confirmation from the original platform documentation. Therefore, `logFC` is reported without interpreting it as a verified log2 fold change.

## 8. Generated Outputs

The analysis generates the following output files.

### Differential Expression Results

- `all_differential_expression.csv` – Complete differential expression analysis results.
- `significant_probes.csv` – Probes satisfying the selected significance criteria.
- `upregulated_probes.csv` – Significant probes with positive logFC values.
- `downregulated_probes.csv` – Significant probes with negative logFC values.
- `annotated_results.csv` – Differential expression results with available platform annotation information.

### Expression and Sample Data

- `clean_expression_matrix.csv` – Preprocessed expression matrix.
- `sample_metadata.csv` – Sample information and experimental conditions.

### Visualizations

- Volcano plot
- Expression heatmap
- MDS plot
- Other exploratory data analysis visualizations

### Research Report

- `final_report.Rmd` – R Markdown source file.
- `final_report.html` – Rendered HTML research report.

## 9. Reproducibility

To reproduce the analysis:

1. Install R and RStudio.
2. Clone or download this repository.
3. Open the `Rice_Drought_ICA3.Rproj` project in RStudio.
4. Install the required R packages if they are not already installed.
5. Run the analysis script located in the `code/` directory.

The dataset is retrieved from NCBI GEO using the GEOquery package.

The analysis results, visualizations, and final report are generated through the R analysis workflow.

## 10. Author

**Undergraduate Research Project**

BSc (Hons) in Information Technology  
Faculty of Applied Science  
University of Vavuniya, Sri Lanka
