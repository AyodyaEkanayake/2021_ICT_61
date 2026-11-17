# Uncovering the Molecular Response of Rice to Drought Stress Using Differential Gene Expression Analysis

## 1. Title and Abstract

### Abstract

Drought stress is an important environmental factor that affects rice growth and productivity. Understanding the molecular responses of rice to drought conditions can help identify genes and biological processes associated with stress responses. This project investigates differential gene expression between control and drought-stressed rice (*Oryza sativa*) samples using publicly available microarray data from the NCBI Gene Expression Omnibus (GEO), dataset GSE115826. The analysis involves data preprocessing, exploratory data analysis, sample relationship assessment, and differential expression analysis using the limma package in R. Statistically significant probes are identified using adjusted P-values and an absolute logFC threshold. Visualizations and annotation information are also used to support the interpretation of the results.

## 2. Introduction

Drought is one of the major environmental stresses affecting agricultural production worldwide. Rice (*Oryza sativa*) is an important food crop, and its growth and productivity can be significantly affected by water scarcity.

Understanding how rice responds to drought stress at the molecular level can provide useful insights into the genes and biological processes involved in stress adaptation.

Differential gene expression analysis is a computational approach used to identify genes or microarray probes that exhibit different expression patterns between experimental conditions.

### Research Aim

The main aim of this project is to investigate the molecular response of rice to drought stress by identifying differentially expressed microarray probes between control and drought-stressed samples.

### Research Objectives

- Obtain publicly available rice gene expression data from NCBI GEO.
- Preprocess and assess the quality of the expression data.
- Explore expression patterns and relationships among samples.
- Identify differentially expressed probes between control and drought-stressed conditions.
- Annotate significant probes using available platform information.
- Visualize and interpret the results to understand potential drought-related expression patterns.

### Research Hypothesis

There are measurable differences in gene expression between control and drought-stressed rice samples.

## 3. Data

### 3.1 Dataset Source

The dataset used in this project was obtained from the National Center for Biotechnology Information (NCBI) Gene Expression Omnibus (GEO).

| Attribute | Description |
|---|---|
| Dataset ID | GSE115826 |
| Source | NCBI Gene Expression Omnibus (GEO) |
| Organism | *Oryza sativa* (Rice) |
| Platform | GPL6864 |
| Data Type | Microarray gene expression data |
| Experimental Conditions | Control and Drought Stress |
| Measurement Type | Agilent Feature Extraction LogRatio values |

**Dataset Link:**  
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115826

### 3.2 Dataset Description

The dataset contains microarray expression measurements from rice samples under control and drought-stressed experimental conditions.

The analysis uses the expression measurements and associated sample metadata obtained from the GEO Series Matrix file.

The main variables used in the analysis include:

- **Probe identifiers:** Unique identifiers representing microarray probes.
- **Expression values:** Original expression measurements obtained from the GEO dataset.
- **Sample identifiers:** Identifiers used to distinguish individual samples.
- **Experimental condition:** Classification of samples as control or drought stress.
- **Annotation information:** Available platform information associated with microarray probes.

### 3.3 Ethical and Licensing Considerations

This project uses publicly available secondary data obtained from NCBI GEO. No new biological samples or personal information were collected for this analysis.

The original dataset source is acknowledged, and the data are used for academic and educational research purposes. Users should refer to the original GEO record and associated data-use conditions when reusing the dataset.

## 4. Methods

### 4.1 Tools and Technologies

The analysis was conducted using R and RStudio.

| Tool / Package | Purpose |
|---|---|
| R | Statistical computing and data analysis |
| RStudio | Development environment for R |
| GEOquery | Retrieving gene expression data from NCBI GEO |
| Biobase | Managing and processing biological expression data |
| limma | Differential expression analysis |
| ggplot2 | Creating statistical visualizations |
| pheatmap | Generating expression heatmaps |
| knitr | Generating report tables and output |
| rmarkdown | Creating the final HTML report |

### 4.2 Analytical Workflow

The project follows the workflow below:

1. Retrieve the GSE115826 dataset from NCBI GEO.
2. Extract expression measurements and sample metadata.
3. Identify control and drought-stressed samples.
4. Perform data quality checks and preprocessing.
5. Conduct exploratory data analysis (EDA).
6. Examine sample relationships using Multidimensional Scaling (MDS).
7. Perform differential expression analysis using limma.
8. Apply empirical Bayes moderation.
9. Correct for multiple hypothesis testing using the Benjamini-Hochberg method.
10. Identify significant, upregulated, and downregulated probes.
11. Retrieve available platform annotation information.
12. Generate visualizations and export analysis results.
13. Prepare the final research report.

### 4.3 Data Preprocessing

The expression values were obtained from the NCBI GEO Series Matrix file.

The dataset contains Agilent Feature Extraction LogRatio measurements. The original GEO values were retained during preprocessing.

Probes containing non-finite values and probes with zero variance across the samples were removed before differential expression analysis.

No additional logarithmic scale conversion was applied. The exact logarithm base of the deposited LogRatio values is not assumed without confirmation from the original platform documentation.

### 4.4 Exploratory Data Analysis

Exploratory Data Analysis (EDA) was performed to understand the characteristics and distribution of the expression data.

The analysis includes:

- Summary statistics of expression measurements.
- Expression distribution comparisons across samples.
- Boxplots for examining sample-level expression distributions.
- Sample distribution visualizations.
- Multidimensional Scaling (MDS) to explore similarities and differences between samples.

These methods help assess the dataset before conducting differential expression analysis.

### 4.5 Differential Expression Analysis

Differential expression analysis was performed using the `limma` package.

A linear modelling approach was used to compare control and drought-stressed samples. Empirical Bayes moderation was then applied to improve the stability of variance estimates.

Multiple hypothesis testing was addressed using the Benjamini-Hochberg correction method.

Probes were classified as significant using the following criteria:

- **Adjusted P-value < 0.05**
- **Absolute logFC (|logFC|) ≥ 1**

Positive logFC values indicate higher expression measurements in the drought-stressed condition relative to the control condition, while negative values indicate lower measurements, according to the model contrast.

The original expression measurements were retained. Therefore, logFC is reported without assuming that it represents a verified log2 fold change.

### 4.6 Probe Annotation

Available platform annotation information was used to provide additional biological context for the identified probes.

Annotation information was combined with the differential expression results where available.

Keyword-based examination of annotation descriptions may be used to explore possible functional themes. Such keyword counts are exploratory observations and do not represent formal functional enrichment analysis.

## 5. Results

The analysis generates statistical results, tables, and visualizations to examine expression differences between control and drought-stressed rice samples.

### 5.1 Differential Expression Results

The following result files are generated:

| File | Description |
|---|---|
| `all_differential_expression.csv` | Complete differential expression analysis results |
| `significant_probes.csv` | Probes satisfying the selected significance criteria |
| `upregulated_probes.csv` | Significant probes with positive logFC values |
| `downregulated_probes.csv` | Significant probes with negative logFC values |
| `annotated_results.csv` | Differential expression results combined with available annotation information |

### 5.2 Expression and Sample Analysis

The exploratory analysis produces visualizations to examine expression distributions and sample relationships.

### 5.3 Generated Figures

The following visualizations are included in the project:

| Figure | Purpose |
|---|---|
| Expression Boxplot | Compares expression distributions across samples |
| Sample Distribution Plot | Shows the distribution of samples across experimental conditions |
| MDS Plot | Explores similarities and differences between samples |
| Volcano Plot | Displays differential expression results using logFC and adjusted significance information |
| Expression Heatmap | Visualizes expression patterns across selected probes and samples |

The figures are saved in the `figures/` directory.

### 5.4 Statistical Interpretation

The differential expression results are interpreted using adjusted P-values and logFC values.

The analysis distinguishes between statistically significant probes with positive and negative logFC values. Available probe annotations are used to support the interpretation of the findings.

The results should be interpreted in the context of the original measurement scale, dataset characteristics, and limitations of microarray-based analysis.

## 6. Discussion

This project applies bioinformatics and statistical analysis techniques to investigate the expression response of rice under drought stress.

The combination of exploratory data analysis, sample relationship assessment, and differential expression analysis provides a structured approach to examining differences between the two experimental conditions.

The limma package provides a suitable statistical framework for microarray differential expression analysis. Empirical Bayes moderation and multiple-testing correction help improve the reliability of the statistical comparisons.

### Limitations

The following limitations should be considered:

- The analysis depends on the quality and completeness of the publicly available GEO dataset.
- Microarray measurements are limited to probes represented on the selected platform.
- Probe annotation information may be incomplete or unavailable for some probes.
- The exact logarithm base of the deposited LogRatio measurements has not been independently confirmed.
- Keyword-based annotation exploration does not establish biological function or statistical enrichment.
- Differential expression results alone do not demonstrate a causal relationship between gene expression and drought tolerance.

The findings should therefore be interpreted as computational evidence of expression differences rather than direct experimental confirmation of biological mechanisms.

### Real-World Implications

Understanding molecular responses to drought stress may contribute to future research on crop stress adaptation, agricultural biotechnology, and the development of drought-resilient rice varieties.

The identified expression patterns can provide a starting point for further biological investigation and experimental validation.

## 7. Conclusion and Future Work

This project investigates the molecular response of rice (*Oryza sativa*) to drought stress using differential gene expression analysis of publicly available microarray data.

The workflow combines data retrieval, preprocessing, exploratory analysis, MDS, differential expression analysis using limma, probe annotation, and visualization.

The analysis provides a computational framework for identifying probes that exhibit expression differences between control and drought-stressed conditions.

### Future Work

Possible improvements and extensions include:

- Performing formal Gene Ontology (GO) and pathway enrichment analysis.
- Investigating the biological functions of significant genes using additional annotation resources.
- Comparing the findings with independent rice drought-stress datasets.
- Validating selected candidate genes using independent experimental evidence.
- Exploring additional statistical and machine learning approaches for identifying drought-response patterns.

## 8. Reproducibility

### 8.1 Project Structure

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

### 8.2 Software Requirements

The project was developed using:

- R 4.6.2
- RStudio 2026.10.1
- Required R packages listed in Section 4.1.

Package versions may differ depending on the environment in which the analysis is executed.

### 8.3 Required R Packages

Install the required packages using the following commands:

```r
install.packages(c(
  "ggplot2",
  "pheatmap",
  "knitr",
  "rmarkdown"
))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c(
  "GEOquery",
  "Biobase",
  "limma"
))
```

### 8.4 Running the Analysis

1. Clone or download this GitHub repository.
2. Open `Rice_Drought_ICA3.Rproj` in RStudio.
3. Ensure the required R packages are installed.
4. Open the `code/analysis.R` script.
5. Run the script to retrieve/process the dataset and generate the analysis outputs.
6. Check the generated files in the `results/` and `figures/` directories.

The dataset is retrieved from NCBI GEO using the GEOquery package.

### 8.5 Generating the Final Report

The final report is written using R Markdown.

Run the following commands in the RStudio Console from the project root directory:

```r
setwd("D:/RLanguage/Rice_Drought_ICA3")

rmarkdown::render(
  "reports/final_report.Rmd",
  output_file = "final_report.html",
  output_dir = "reports"
)
```

The rendered HTML report will be saved in the `reports/` directory.

## 9. References

1. National Center for Biotechnology Information (NCBI). Gene Expression Omnibus (GEO): GSE115826. Available at: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115826

2. R Core Team. R: A Language and Environment for Statistical Computing. R Foundation for Statistical Computing, Vienna, Austria. Available at: https://www.r-project.org/

3. RStudio Team. RStudio: Integrated Development Environment for R. Posit Software. Available at: https://posit.co/products/open-source/rstudio/

4. Ritchie, M. E., Phipson, B., Wu, D., Hu, Y., Law, C. W., Shi, W., and Smyth, G. K. (2015). limma powers differential expression analyses for RNA-sequencing and microarray studies. *Nucleic Acids Research*, 43(7), e47. https://doi.org/10.1093/nar/gkv007

5. Davis, S., and Meltzer, P. S. (2007). GEOquery: a bridge between the Gene Expression Omnibus (GEO) and BioConductor. *Bioinformatics*, 23(14), 1846–1847. https://doi.org/10.1093/bioinformatics/btm254

6. Wickham, H. ggplot2: Elegant Graphics for Data Analysis. Springer. https://ggplot2.tidyverse.org/

7. Kolde, R. pheatmap: Pretty Heatmaps. R package. https://cran.r-project.org/package=pheatmap

---

**Author:** Undergraduate Research Project  
**Degree:** BSc (Hons) in Information Technology  
**Faculty:** Faculty of Applied Science  
**University:** University of Vavuniya, Sri Lanka