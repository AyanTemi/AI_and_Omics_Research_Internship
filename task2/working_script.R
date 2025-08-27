# Define a function to classify genes as unregulated, downregulated or not significant
classify_gene <- function(logFC, padj) {
  if (logFC > 1 & padj < 0.05) {
    return("Upregulated")
  } else if (logFC < -1 & padj < 0.05) {
    return("Downregulated")
  } else {
    return("Not_Significant")
  }
}

# Define input and output folders
input_dir <- "Raw_Data" 
output_dir <- "Results"

#create output folder
if(!dir.exists(output_dir)){
  dir.create(output_dir)
}

# List which files to process
files_to_process <- c("DEGs_Data_1.csv", "DEGs_Data_2.csv") 


# For each file within a loop:
#       - Apply the function to classify genes in a for-loop to process both datasets
#       - Replace missing padj values with 1
#       - Add a new column 'status'
#       - save processed files into Results folder
#       - Print summary counts of significant, upregulated, and downregulated genes
#       - Use table() for summaries

for (file_names in files_to_process) {
  cat("\nProcessing:", file_names, "\n")
  
  input_file_path <- file.path(input_dir, file_names)
  
  # Import dataset
  data <- read.csv(input_file_path, header = TRUE)
  cat("File imported. Checking for missing values...\n")
  
  # handling missing values
  data$padj[is.na(data$padj)] <- 1
  
  # Create an empty column for status
  data$status <- NA
  
  # Apply classification (vectorised with mapply)
  data$status <- mapply(classify_gene, data$logFC, data$padj)
  
  # Save processed file
  out_file <- paste0("Results/processed_", file_names)
  write.csv(data, out_file, row.names = FALSE)
  
  # Print summary
  cat("\nSummary for:", file_names, "\n")
  print(table(data$status))
}

#Save the entire R workspace
save.image(file = "TemitopeAyano_Class_2_Assignment.RData")
