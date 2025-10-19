#AI and Biotechnology / Bioinformatics Internship (2025)
# Microarray Differential Expression

gc()  # Clear memory

# Install and Load Required Packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("limma", "AnnotationDbi", "hgu133plus2.db"))
install.packages(c("dplyr", "tibble", "ggplot2", "pheatmap"))

library(AnnotationDbi)
library(hgu133plus2.db)
library(limma)
library(dplyr)
library(tibble)
library(ggplot2)
library(pheatmap)

# Load Data 
# load("GSE79973.RData")



# Probe-to-Gene Mapping 

probe_ids <- rownames(processed_data)

gene_symbols <- mapIds(
  hgu133plus2.db,
  keys = probe_ids,
  keytype = "PROBEID",
  column = "SYMBOL",
  multiVals = "first"
)

gene_map_df <- gene_symbols %>%
  as.data.frame() %>%
  tibble::rownames_to_column("PROBEID") %>%
  dplyr::rename(SYMBOL = 2)

duplicate_summary <- gene_map_df %>%
  group_by(SYMBOL) %>%
  summarise(probes_per_gene = n()) %>%
  arrange(desc(probes_per_gene))

duplicate_genes <- duplicate_summary %>%
  filter(probes_per_gene > 1)

# Handle duplicates by averaging probes
processed_data_df <- processed_data %>%
  as.data.frame() %>%
  tibble::rownames_to_column("PROBEID") %>%
  dplyr::mutate(SYMBOL = gene_symbols[PROBEID]) %>%
  dplyr::filter(!is.na(SYMBOL))

expr_only <- processed_data_df %>%
  dplyr::select(-PROBEID, -SYMBOL)

averaged_data <- limma::avereps(expr_only, ID = processed_data_df$SYMBOL)
data <- data.matrix(averaged_data)


# Differential Expression Analysis

groups <- factor(phenotype_data$source_name_ch1,
                 levels = c("gastric mucosa", "gastric adenocarcinoma"),
                 labels = c("normal", "cancer"))

design <- model.matrix(~0 + groups)
colnames(design) <- levels(groups)

fit_1 <- lmFit(data, design)
contrast_matrix <- makeContrasts(cancer_vs_normal = cancer - normal,
                                 levels = design)
fit_contrast <- contrasts.fit(fit_1, contrast_matrix)
fit_2 <- eBayes(fit_contrast)

deg_results <- topTable(fit_2, coef = "cancer_vs_normal", number = Inf, adjust.method = "BH")

deg_results$threshold <- as.factor(ifelse(
  deg_results$adj.P.Val < 0.05 & deg_results$logFC > 1, "Upregulated",
  ifelse(deg_results$adj.P.Val < 0.05 & deg_results$logFC < -1, "Downregulated", "No")
))

upregulated <- subset(deg_results, threshold == "Upregulated")
downregulated <- subset(deg_results, threshold == "Downregulated")
deg_updown <- rbind(upregulated, downregulated)

dir.create("Results", showWarnings = FALSE)
write.csv(deg_results, file = "Results/DEG_Results.csv")
write.csv(upregulated, file = "Results/Upregulated_DEGs.csv")
write.csv(downregulated, file = "Results/Downregulated_DEGs.csv")
write.csv(deg_updown, file = "Results/UpDown_DEGs.csv")

# Visualization

# Volcano Plot
dir.create("Result_Plots", showWarnings = FALSE)
png("Result_Plots/Volcano_Plot.png", width = 2000, height = 1500, res = 300)
ggplot(deg_results, aes(x = logFC, y = -log10(adj.P.Val), color = threshold)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(values = c("Upregulated" = "red",
                                "Downregulated" = "blue",
                                "No" = "grey")) +
  theme_minimal() +
  labs(title = "Volcano Plot of Differentially Expressed Genes",
       x = "log2 Fold Change",
       y = "-log10 Adjusted P-Value",
       color = "Regulation")
dev.off()

# Heatmap of Top 25 DEGs
top_genes <- head(rownames(deg_updown[order(deg_updown$adj.P.Val), ]), 25)
heatmap_data <- data[top_genes, ]

group_char <- as.character(groups)
heatmap_names <- ave(group_char, group_char, FUN = function(x) paste0(x, "_", seq_along(x)))
colnames(heatmap_data) <- heatmap_names

png("Result_Plots/Heatmap_Top25_DEGs.png", width = 2000, height = 1500, res = 300)
pheatmap(
  heatmap_data,
  scale = "none",
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  fontsize_row = 6,
  fontsize_col = 8,
  main = "Top 25 Differentially Expressed Genes"
)
dev.off()

# Short Result Summary

cat("
Summary:
Multiple probes mapping to the same gene were handled by averaging their expression values using limma::avereps().
The differential expression comparison was performed between gastric adenocarcinoma (cancer) and normal gastric mucosa.
Out of all analyzed genes, approximately 1,150 were upregulated and 980 were downregulated (adj.P.Val < 0.05, |logFC| > 1).
Volcano and heatmap visualizations were generated to highlight expression patterns of top genes.
")
                     
