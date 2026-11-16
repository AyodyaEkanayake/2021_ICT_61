# 1. PROJECT SETUP
# If the script is opened from inside the code/ folder, move to the project root
if (basename(getwd()) %in% c("code", "reports")) setwd("..")

# Create project folders
folders <- c(
  "data",
  "results",
  "figures",
  "reports"
)

for (folder in folders) {
  dir.create(
    folder,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

# Set package library location
package_path <- "D:/RLanguage/R_Packages"

# Create package folder if needed
dir.create(
  package_path,
  recursive = TRUE,
  showWarnings = FALSE
)

# Set library path
.libPaths(c(package_path, .libPaths()))

# Check library path
.libPaths()


# 2. PACKAGES
# Packages are only LOADED here (not re-installed on every run).
# If a package is missing, run the install commands below ONCE in the Console:
#
#   install.packages("BiocManager", lib = package_path)
#   BiocManager::install(c("GEOquery", "limma", "Biobase"), lib = package_path)
#   install.packages(c("ggplot2", "pheatmap"), lib = package_path)

library(GEOquery)
library(limma)
library(Biobase)
library(ggplot2)
library(pheatmap)

# Display package versions
packageVersion("GEOquery")
packageVersion("limma")
packageVersion("Biobase")
packageVersion("pheatmap")


# 3. LOAD THE GEO DATASET
# Set dataset ID
dataset_id <- "GSE115826"

# Use the Series Matrix file downloaded from NCBI if it is available
# (the online getGEO() download can fail on some networks).
# https://ftp.ncbi.nlm.nih.gov/geo/series/GSE115nnn/GSE115826/matrix/GSE115826_series_matrix.txt.gz
local_matrix <- "data/GSE115826_series_matrix.txt.gz"

if (file.exists(local_matrix)) {
  
  gse <- getGEO(
    filename = local_matrix,
    GSEMatrix = TRUE,
    getGPL = FALSE
  )
  
} else {
  
  # Download dataset from NCBI GEO
  gse <- getGEO(
    dataset_id,
    GSEMatrix = TRUE,
    getGPL = FALSE,
    destdir = "data"
  )
  
}

# The online download returns a LIST of platforms,
# the local file returns a single ExpressionSet. Handle both.
if (is.list(gse)) {
  eset <- gse[[1]]
} else {
  eset <- gse
}

# Display dataset information
eset


# 4. EXTRACT EXPRESSION DATA, SAMPLE INFORMATION AND PROBE INFORMATION

# Extract gene expression values
expression_data <- exprs(eset)

# Extract sample information
sample_data <- pData(eset)

# Extract probe information
probe_data <- fData(eset)

# Display first 6 expression rows
head(expression_data)

# Check expression matrix dimensions
dim(expression_data)

# Check sample names
colnames(expression_data)

# Check sample information
head(sample_data)

# Check probe information
head(probe_data)

# Save the downloaded dataset
saveRDS(
  eset,
  "data/GSE115826_dataset.rds"
)

# Save expression data
write.csv(
  expression_data,
  "data/expression_data.csv"
)

# Save sample information
write.csv(
  sample_data,
  "data/sample_information.csv"
)

# Check sample metadata
sample_data$title
sample_data$geo_accession
sample_data$characteristics_ch1


# 5. IDENTIFY THE CONTROL AND DROUGHT GROUPS
# Extract sample titles
sample_titles <- sample_data$title

# Identify sample conditions
condition <- ifelse(
  grepl("drought", sample_titles, ignore.case = TRUE),
  "Drought",
  "Control"
)

# Convert into factor
condition <- factor(
  condition,
  levels = c("Control", "Drought")
)

# Display groups
print(condition)

# Count samples in each group (expected: 4 Control, 4 Drought)
table(condition)
stopifnot(all(table(condition) == 4))

# Add condition to the sample metadata
sample_data$condition <- condition

# Display updated metadata
sample_data[, c("title", "condition")]


# 6. DATA QUALITY CHECKS AND PREPROCESSING
# Missing values (total and per sample)
sum(is.na(expression_data))
colSums(is.na(expression_data))

# Duplicate probe IDs
sum(duplicated(rownames(expression_data)))

# Non-finite (NA / NaN / Inf) values
sum(!is.finite(expression_data))

# Summary of expression values
summary(as.vector(expression_data))

# Remove probes with missing or infinite values
valid_rows <- apply(
  expression_data,
  1,
  function(x) all(is.finite(x))
)

# Remove zero-variance probes
variable_rows <- apply(
  expression_data,
  1,
  var
) > 0

# Keep valid probes
clean_expression <- expression_data[
  valid_rows & variable_rows,
]

# Convert log10 ratios into log2 ratios.
# The GEO Series Matrix values are the Agilent Feature Extraction "LogRatio"
# column, which is log10(Cy5 / Cy3). Multiplying by log2(10) converts them to
# log2 ratios, so that logFC is on the log2 scale (log2(x) = log10(x) * log2(10)).
clean_expression <- clean_expression * log2(10)

# Check dimensions (probes x samples)
dim(clean_expression)


# 7. EXPLORATORY DATA ANALYSIS (EDA)
# ---- 7.1 Sample distribution ----

# Count each group
group_counts <- data.frame(
  Group = names(table(condition)),
  Samples = as.numeric(table(condition))
)

# Create bar chart
sample_plot <- ggplot(
  group_counts,
  aes(x = Group, y = Samples, fill = Group)
) +
  geom_bar(
    stat = "identity"
  ) +
  scale_fill_manual(
    values = c("Control" = "#5B7FA6", "Drought" = "#C0504D")
  ) +
  labs(
    title = "Distribution of Rice Samples",
    x = "Experimental Group",
    y = "Number of Samples"
  ) +
  theme_minimal()

print(sample_plot)

# Save sample distribution plot
ggsave(
  "figures/sample_distribution.png",
  plot = sample_plot,
  width = 7,
  height = 5,
  dpi = 300
)


# ---- 7.2 Boxplot - sample expression distribution ----

# Convert expression matrix into a data frame
expression_df <- as.data.frame(clean_expression)

# Boxplot colours follow the experimental group
box_colours <- c("#5B7FA6", "#C0504D")[as.numeric(condition)]

# Draw boxplot on screen
boxplot(
  expression_df,
  main = "Expression Distribution Across Samples",
  xlab = "Samples",
  ylab = "Log2 Expression Ratio",
  col = box_colours,
  las = 2
)

# Save boxplot
png(
  "figures/expression_boxplot.png",
  width = 2000,
  height = 1500,
  res = 300
)

boxplot(
  expression_df,
  main = "Expression Distribution Across Samples",
  xlab = "Samples",
  ylab = "Log2 Expression Ratio",
  col = box_colours,
  las = 2
)

dev.off()


# ---- 7.3 PCA / MDS analysis ----
# Check similarities between microarray samples

# Draw MDS plot on screen
plotMDS(
  clean_expression,
  labels = condition,
  col = as.numeric(condition),
  main = "MDS Plot of Rice Samples"
)

# Save MDS plot
png(
  "figures/mds_plot.png",
  width = 2000,
  height = 1500,
  res = 300
)

plotMDS(
  clean_expression,
  labels = condition,
  col = as.numeric(condition),
  main = "MDS Plot of Rice Samples"
)

dev.off()


# 8. DIFFERENTIAL GENE EXPRESSION ANALYSIS (limma)
# Finds: upregulated probes, downregulated probes, non-significant probes

# ---- 8.1 Experimental design ----

# Create design matrix
design <- model.matrix(
  ~ 0 + condition
)

# Rename columns
colnames(design) <- c(
  "Control",
  "Drought"
)

# Display design
design

# ---- 8.2 Contrast: Drought - Control ----
# Positive logFC = higher expression under drought
# Negative logFC = lower expression under drought

contrast_matrix <- makeContrasts(
  Drought - Control,
  levels = design
)

# Display contrast
contrast_matrix

# ---- 8.3 Linear model + empirical Bayes ----

# Fit linear model
fit <- lmFit(
  clean_expression,
  design
)

# Apply contrast
fit2 <- contrasts.fit(
  fit,
  contrast_matrix
)

# Apply empirical Bayes
fit2 <- eBayes(fit2)

# ---- 8.4 Extract all results ----
results <- topTable(
  fit2,
  coef = 1,
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)

# Display first 10 results
head(results, 10)


# 9. SIGNIFICANT, UPREGULATED AND DOWNREGULATED PROBES
# Criteria: adjusted P-value < 0.05 AND |log2 fold change| >= 1

significant <- results[
  results$adj.P.Val < 0.05 &
    abs(results$logFC) >= 1,
]

# Display significant results
head(significant)

# Count significant probes
nrow(significant)

# Upregulated probes
upregulated <- significant[
  significant$logFC > 0,
]

# Downregulated probes
downregulated <- significant[
  significant$logFC < 0,
]

# Display counts
nrow(upregulated)
nrow(downregulated)


# 10. VOLCANO PLOT
# Prepare results for volcano plot
results$ProbeID <- rownames(results)

results$Status <- "Not Significant"

results$Status[
  results$adj.P.Val < 0.05 &
    results$logFC >= 1
] <- "Upregulated"

results$Status[
  results$adj.P.Val < 0.05 &
    results$logFC <= -1
] <- "Downregulated"

# Avoid infinite values
results$MinusLogP <- -log10(
  pmax(results$adj.P.Val, .Machine$double.xmin)
)

# Create volcano plot
volcano <- ggplot(
  results,
  aes(
    x = logFC,
    y = MinusLogP,
    color = Status
  )
) +
  geom_point(alpha = 0.7, size = 1.5) +
  scale_color_manual(
    values = c(
      "Upregulated"     = "#C0504D",
      "Downregulated"   = "#2F5D8A",
      "Not Significant" = "#A6AEB8"
    )
  ) +
  geom_vline(
    xintercept = c(-1, 1),
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed"
  ) +
  labs(
    title = "Volcano Plot of Rice Drought Stress",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value"
  ) +
  theme_minimal()

print(volcano)

# Save the volcano plot
ggsave(
  "figures/volcano_plot.png",
  plot = volcano,
  width = 8,
  height = 6
)


# 11. HEATMAP (top 50 most variable probes)
# Calculate probe variance FIRST (it is needed in both branches below)
probe_variance <- apply(
  clean_expression,
  1,
  var
)

# Select probes for heatmap (fewer than 50 probes are handled safely)
n_heatmap <- min(50, length(probe_variance))

selected_probes <- names(
  sort(probe_variance, decreasing = TRUE)
)[1:n_heatmap]

# Prepare sample annotation
sample_annotation <- data.frame(
  Condition = condition
)

rownames(sample_annotation) <- colnames(clean_expression)

# Draw heatmap
pheatmap(
  clean_expression[selected_probes, ],
  scale = "row",
  annotation_col = sample_annotation,
  show_rownames = FALSE,
  main = "Rice Drought Stress Gene Expression Heatmap"
)

# Save heatmap as PNG
pheatmap(
  clean_expression[selected_probes, ],
  scale = "row",
  annotation_col = sample_annotation,
  show_rownames = FALSE,
  main = "Rice Drought Stress Gene Expression Heatmap",
  filename = "figures/expression_heatmap.png",
  width = 10,
  height = 8
)

#----------------------------------------------------------------------------------------
# Download platform information (uses data/GPL6864.soft.gz if already present)
platform <- getGEO(
  "GPL6864",
  destdir = "data"
)

# Extract annotation table
platform_data <- Table(platform)

# Display annotation columns
names(platform_data)

# Display first 6 rows
head(platform_data)

# Add probe IDs
significant$ProbeID <- rownames(significant)

# Match probe IDs with platform annotation
annotated_results <- merge(
  significant,
  platform_data,
  by.x = "ProbeID",
  by.y = "ID",
  all.x = TRUE
)

# Display annotated results
head(annotated_results)


# 13. SAVE ALL RESULTS
# Save all differential expression results
write.csv(
  results,
  "results/all_differential_expression.csv",
  row.names = FALSE
)

# Save significant probes
write.csv(
  significant,
  "results/significant_probes.csv",
  row.names = FALSE
)

# Save upregulated probes
write.csv(
  upregulated,
  "results/upregulated_probes.csv"
)

# Save downregulated probes
write.csv(
  downregulated,
  "results/downregulated_probes.csv"
)

# Save annotated results
write.csv(
  annotated_results,
  "results/annotated_results.csv",
  row.names = FALSE
)

# Save clean expression matrix
write.csv(
  clean_expression,
  "results/clean_expression_matrix.csv"
)

# Save sample metadata
write.csv(
  sample_data,
  "results/sample_metadata.csv"
)


# 14. FINAL RESULTS

cat("Total samples:", ncol(clean_expression), "\n")

cat("Total probes:", nrow(clean_expression), "\n")

cat("Significant probes:", nrow(significant), "\n")

cat("Upregulated probes:", nrow(upregulated), "\n")

cat("Downregulated probes:", nrow(downregulated), "\n")
