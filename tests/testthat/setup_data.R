library(here) #FIXME remove this!!!
## FIXME This needs to be data(example_seurat) for example.
seurat_object <- readRDS(here::here("data","example_data.Rds")) # Load data
sample_ids <- seurat_object[[]]$sample # Get sample IDs for each cell

# get full subset
full_pca_embeddings <- seurat_object@reductions$pca@cell.embeddings # Extract PCA embeddings from Seurat object
full_pca_embeddings_subset <- full_pca_embeddings[,1:10] # Warning: Never subsample scVI embeddings

# subset one sample to 49 cells to check cell count-based warnings
subset_sample_id <- sample_ids[20]
sample_drop_indices <- sample(which(seurat_object[[]]$sample==subset_sample_id),(500-49),replace=FALSE)
reliability_pca_embeddings_subset <- full_pca_embeddings_subset[-sample_drop_indices,]
reliability_sample_ids <- sample_ids[-sample_drop_indices]

### Keep this for compatibility with old tests, at least for now
data(example_data)
sub_pat <- sample(unique(example_data$patient_id),3)
sub_data <- example_data[example_data$patient_id %in% sub_pat,]

