library(here)
seurat_object <- readRDS(here::here("data","testing_data.Rds")) # Load data             
sample_ids <- seurat_object[[]]$sample # Get sample IDs for each cell     
pca_embeddings <- seurat_object@reductions$pca@cell.embeddings # Extract PCA embeddings from Seurat object
pca_embeddings_subset <- pca_embeddings[,1:10] # Warning: Never subsample scVI embeddings
