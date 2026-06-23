#This is how to load and process data in preparation for dWUI calculation
library(dplyr)
library(Seurat)
library(patchwork)

# Path that contains your GSM folders
data_path <- "~/amramujc/Paper1_GALNT2/Sample_Folders"

# Get only folders that start with GSM
sample_dirs <- list.dirs(data_path, recursive = FALSE, full.names = TRUE)
sample_dirs <- sample_dirs[grepl("^GSM", basename(sample_dirs))]

# Verify which folders will be processed
cat("Found these 10X data directories:\n")
print(basename(sample_dirs))

processed_list <- list()

for (dir in sample_dirs) {
  sample_id <- basename(dir)
  sample_id <- sub("_.*", "", basename(dir))
  cat("\nProcessing:", sample_id, "\n")
  sample_id <- sub("_.*", "", basename(dir))
  
  # Read 10X data
  data <- Read10X(data.dir = dir)
  seurat <- CreateSeuratObject(
    counts = data,
    project = sample_id,
    min.cells = 3,
    min.features = 200
  )
  
  # Calculate %MT
  seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")
  
  # QC
  seurat <- subset(seurat,
                   subset = nFeature_RNA > 200 &
                     nFeature_RNA < 8000 &
                     percent.mt < 10)
  
  # Add sample metadata
  seurat$Sample <- sample_id
  
  # Normalize,  Scaling
  seurat <- NormalizeData(seurat, normalization.method = "LogNormalize", scale.factor = 10000)
  seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000)
  
  seurat <- ScaleData(seurat, features = rownames(seurat))
  seurat <- RunPCA(seurat, features = VariableFeatures(seurat))
  
  # Neighbors / Clusters / UMAP (optional pre-merge)
  seurat <- FindNeighbors(seurat, dims = 1:10)
  seurat <- FindClusters(seurat, resolution = 0.5)
  seurat <- RunUMAP(seurat, dims = 1:10)
  
  # Save RDS
  saveRDS(seurat, file = paste0(sample_id, "_processed.rds"))
  processed_list[[sample_id]] <- seurat
}

cat("\nComplete: processed and saved all samples!\n")
