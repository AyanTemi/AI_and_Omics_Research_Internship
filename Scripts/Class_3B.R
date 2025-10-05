# ---------------------------------------------------------------------
#                   Microarray Data Analysis Assignment
# =====================================================================

# Install and Load Required Packages 

# Install Bioconductor (provides R packages for analyzing omics data)

if (!requireNamespace("BiocManager", quietly = TRUE)) 
  install.packages("BiocManager")

# Install Bioconductor packages
BiocManager::install(c("GEOquery","affy","arrayQualityMetrics"))

# Install CRAN packages for data manipulation
install.packages("dplyr")

# Load Required Libraries
library(GEOquery)             # Download GEO datasets (series matrix, raw CEL files)
library(affy)                 # Pre-processing of Affymetrix microarray data (RMA normalization)
library(arrayQualityMetrics)  # QC reports for microarray data
library(dplyr)                # Data manipulation

# -------------------------------------
#### Download Series Matrix Files ####
# -------------------------------------

# Series matrix files are preprocessed text files containing 
# expression values, sample annotations, and probe information.

gse_data <- getGEO("GSE79973", GSEMatrix = TRUE)

# Extract expression data matrix (genes/probes × samples)
# Rows corresponds to probes and columns corresponds to samples
expression_data <- exprs(gse_data$GSE79973_series_matrix.txt.gz)


# Extract feature (probe annotation) data
# Rows corresponds to probes and columns corresponds to samples
feature_data <-  fData(gse_data$GSE79973_series_matrix.txt.gz)


# Extract phenotype (sample metadata) data
# Rows corresponds to samples and columns corresponds to probes
phenotype_data <-  pData(gse_data$GSE79973_series_matrix.txt.gz)

# Check missing values in sample annotation
sum(is.na(phenotype_data$source_name_ch1)) 

# Fetch GEO supplementry files
getGEOSuppFiles("GSE79973", baseDir = "Data", makeDirectory = TRUE)

# Untar CEL files if compressed as .tar
untar("Data/GSE79973/GSE79973_RAW.tar", exdir = "Data/CEL_Files")

# Read CEL files into R as an AffyBatch object
raw_data <- ReadAffy(celfile.path = "Data/CEL_Files")

raw_data   # Displays basic information about the dataset
# annotation= hgu133plus2

# ---------------------------------------------------------------------------
# 1. Perform quality control before and after normalization and 
# check whether any arrays are flagged as outliers. 
# note down how many you found before and after normalization
# ---------------------------------------------------------------------------

arrayQualityMetrics(expressionset = raw_data,
                    outdir = "Results/QC_Raw_Data",
                    force = TRUE,
                    do.logtransform = TRUE)

# Check whether any arrays are flagged as outliers.
# Yes, 5 arrays before normalizaton and 1 after normalization

# ---------------------------------------------------------------------------
# 2. Normalize the data and then apply filtering to remove low-intensity probes 
# and note how many transcripts remain. 
# ---------------------------------------------------------------------------

normalized_data <- rma(raw_data)

# QC after data normalization 
arrayQualityMetrics(expressionset = normalized_data,
                    outdir = "Results/QC_Normalized_Data",
                    force = TRUE)

# Extract normalized expression values into a data frame
processed_data <- as.data.frame(exprs(normalized_data))

dim(processed_data)   # Dimensions: number of probes × number of samples

# Filter Low-Variance Transcripts

# Calculate median intensity per probe across samples
row_median <- rowMedians(as.matrix(processed_data))

# Visualize distribution of probe median intensities
hist(row_median,
     breaks = 100,
     freq = FALSE,
     main = "Median Intensity Distribution")

# Set a threshold to remove low variance probes
threshold <- 3.5 
abline(v = threshold, col = "black", lwd = 2) 

# Select probes above threshold
indx <- row_median > threshold 
filtered_data <- processed_data[indx, ] 

# Rename filtered expression data with sample metadata
colnames(filtered_data) <- rownames(phenotype_data)

# Overwrite processed data with filtered dataset
processed_data <- filtered_data 


# ---------------------------------------------------------------------------
# 3. Use the phenotype information to define your target groups and re-label them 
# (e.g normal vs cancer)
# ---------------------------------------------------------------------------
class(phenotype_data$source_name_ch1) 

# Define experimental groups (normal vs cancer)
groups <- factor(phenotype_data$source_name_ch1,
                 levels = c("gastric mucosa", "gastric adenocarcinoma"),
                 label = c("normal", "cancer"))

class(groups)
levels(groups)
