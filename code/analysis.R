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

################################################
# Insrall heatmap package
package_path <- "D:/RLanguage/R_Packages"

dir.create(
  package_path,
  recursive = TRUE,
  showWarnings = FALSE
)

install.packages(
  "pheatmap",
  lib = package_path,
  repos = "https://cloud.r-project.org",
  type = "binary"
)

# load the packages
.libPaths(c(package_path, .libPaths()))

library(pheatmap)

# verify installation
packageVersion("pheatmap")



library(GEOquery)
library(limma)
library(Biobase)
library(ggplot2)
library(pheatmap)

packageVersion("GEOquery")

packageVersion("limma")

packageVersion("Biobase")


# Load GEOquery
library(GEOquery)

# Set dataset ID
dataset_id <- "GSE115826"

# Download dataset from NCBI GEO
gse <- getGEO(
  dataset_id,
  GSEMatrix = TRUE,
  getGPL = FALSE,
  destdir = "data"
)

# Check number of platforms
length(gse)


###############################################

#get dataset object
# Select the first platform
eset <- gse[[1]]

# Display dataset information
eset


###############################################
# get expression data
# Extract gene expression values
expression_data <- exprs(eset)

# Extract sample information
sample_data <- pData(eset)

# Extract probe information
probe_data <- fData(eset)

# Display first 6 expression rows
head(expression_data)


###########################################
# Check the dataset size
# Check expression matrix dimensions
dim(expression_data)

# Check sample names
colnames(expression_data)

# Check sample information
head(sample_data)

# Check probe information
head(probe_data)


############################################
#Save the dataset
# Save downloaded dataset
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


##############################################
#Check sample metadata
# Display sample titles
sample_data$title

# Display GEO sample IDs
sample_data$geo_accession

# Display sample characteristics
sample_data$characteristics_ch1


#############################################
# Identify the control and drought grps
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




############################################
# Check grp count
# Count samples in each group
table(condition)


#############################################
# Add sample metadata - Final stage of identifying experimental groups
sample_data$condition <- condition

# Display updated metadata
sample_data[, c("title", "condition")]


################################################
#Data Preprocessing - modify before use
# Check the missing values
# Check missing values
sum(is.na(expression_data))

# Check missing values in each sample
colSums(is.na(expression_data))


#################################################
#Check for infinite values
# Check infinite values

###################################################
## Check duplicate probe IDs
sum(duplicated(rownames(expression_data)))
sum(!is.finite(expression_data))

#########################################
# Summary of expression values
summary(as.vector(expression_data))


#################################################
# Data Preprocessing
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

# Convert log10 ratios into log2 ratios
clean_expression <- clean_expression * log2(10)

# Check dimensions
dim(clean_expression)


#######################################################
# Exploratory Data Analysis (EDA)
# Sample Distribution

# Count each group
group_counts <- data.frame(
  Group = names(table(condition)),
  Samples = as.numeric(table(condition))
)

# Create bar chart
ggplot(
  group_counts,
  aes(x = Group, y = Samples, fill = Group)
) +
  geom_bar(
    stat = "identity"
  ) +
  labs(
    title = "Distribution of Rice Samples",
    x = "Experimental Group",
    y = "Number of Samples"
  ) +
  theme_minimal()



##################################################
# Boxplot – Sample Expression Distribution

# Convert expression matrix into a data frame
expression_df <- as.data.frame(clean_expression)

# Create boxplot
boxplot(
  expression_df,
  main = "Expression Distribution Across Samples",
  xlab = "Samples",
  ylab = "Log2 Expression Ratio",
  las = 2
)


##############################################
# PCA / MDS Analysis
# Check similarities between Microarray data
# Create MDS plot
plotMDS(
  clean_expression,
  labels = condition,
  col = as.numeric(condition),
  main = "MDS Plot of Rice Samples"
)


####################################################
#Differential Gene Expression Analysis (Main Analysis)
# Find - Upregulated probes
# Downregulated probes
# Non-significant probes

# make Experimental Design 
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


###############################################
# Compare Control and Drought
# Create contrast
contrast_matrix <- makeContrasts(
  Drought - Control,
  levels = design
)

# Display contrast
contrast_matrix


#####################################################
# # Fit linear model
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


#################################################
# Get significant results 
# Extract all results
results <- topTable(
  fit2,
  coef = 1,
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)

# Display first 10 results
head(results, 10)


###############################################
# Identify significant probes
significant <- results[
  results$adj.P.Val < 0.05 &
    abs(results$logFC) >= 1,
]

# Display significant results
head(significant)

# Count significant probes
nrow(significant)

###################################################
# Upregulated and Downregulated Probes
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


###############################################
# Create volcano plot

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



################################################################
# save the volcano 

ggsave(
  "figures/volcano_plot.png",
  plot = volcano,
  width = 8,
  height = 6
)


################################################################
# Create the heatmap

# Select probes for heatmap
if (nrow(significant) >= 2) {
  
  selected_probes <- head(
    names(sort(probe_variance, decreasing = TRUE)),
    50
  )
  
} else {
  
  probe_variance <- apply(
    clean_expression,
    1,
    var
  )
  
  selected_probes <- names(
    sort(
      probe_variance,
      decreasing = TRUE
    )
  )[1:50]
  
}

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



#################################################################

# save the heatmap
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


###############################################################

# Microarray Probes + Gene Annotations

# Download Platform Annotation 
# Download platform information
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




#####################################################

# Significant Results + Annotation 
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



#####################################################
# Save all the reuslts
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



############################################################
# Final results

cat("Total samples:", ncol(clean_expression), "\n")

cat("Total probes:", nrow(clean_expression), "\n")

cat("Significant probes:", nrow(significant), "\n")

cat("Upregulated probes:", nrow(upregulated), "\n")

cat("Downregulated probes:", nrow(downregulated), "\n")