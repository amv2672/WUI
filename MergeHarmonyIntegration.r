#Merge and Harmony Integration

library(Seurat)
library(dplyr)
library(patchwork)
library(harmony)

#Load
rds_files <- list.files("~/amramujc/Kim", pattern = "_processed\\.rds$", full.names = TRUE)

# Read all RDS files into a named list
sample_list <- lapply(rds_files, readRDS)
# Metadata table
metadata_table <- data.frame(
  Sample = c("GSM5944923", "GSM5944922"),
  Condition = c("Cancer", "Normal"),
  PatientID = c("P01","P01")
)

# Annotate each Seurat object and update barcodes
for (i in seq_along(sample_list)) {
  
  sample <- sample_list[[i]]
  
  # Extract sample name from file name
  sample_id <- gsub("_processed\\.rds$", "", basename(rds_files[i]))
  
  # Find metadata for this sample
  meta <- metadata_table %>% filter(Sample == sample_id)
  
  # Add metadata columns
  sample$Sample <- meta$Sample
  sample$Condition <- meta$Condition
  sample$PatientID <- meta$PatientID
  
  # Build ID for Harmony or barcode suffix (C09/N09)
  patient_num <- gsub("P", "", meta$PatientID)       # remove "P"
  ID <- paste0(substr(meta$Condition, 1, 1), patient_num)  # "C09" or "N09"
  sample$ID <- ID
  
  # Append ID as suffix to cell barcodes
  colnames(sample) <- paste0(colnames(sample), "-", ID)
  
  # Save back
  sample_list[[i]] <- sample
}
# Merge
merged <- merge(sample_list[[1]], y = sample_list[-1], project = "MergedSamples")

# Normalize
merged <- NormalizeData(merged)
merged <- FindVariableFeatures(merged)
merged <- ScaleData(merged)
merged <- RunPCA(merged)
saveRDS(merged, file = "MergedSamples_PriorToIntegration_Kim_P01_NT_TC_1207.rds")
# Harmony integration
# Integrate across "Sample" to remove batch/sample effects
merged <- IntegrateLayers(
  object = merged,
  method = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction = "harmony",
  verbose = FALSE
)

# UMAP
merged <- RunUMAP(merged, reduction = "harmony", dims = 1:30)
merged <- FindNeighbors(merged, reduction = "harmony", dims = 1:30)
merged <- FindClusters(merged, resolution = 0.5)

# Save
saveRDS(merged, file = "MergedSamples_HarmonyIntegrated_Kim_P01_NT_TC_1207.rds")

cat("\nHarmony integration complete! Saved as 'MergedSamples_HarmonyIntegrated_Kim_P01_NT_TC_1207.rds'\n")
