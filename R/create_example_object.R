#library(dplyr)
#library(here)
#library(Seurat)
# set.seed(2)
# print("Seed set")

# Load the full pre-processed 143 sample COVID PBMC data
#seurat_object <- readRDS(here::here("..","GloScope_analysis","data","Processed_Datasets","stephensonCOVIDPBMC","stephensonCOVIDPBMC_default","stephensonCOVIDPBMC_default.Rds"))
# we only keep samples from the Cambridge sequencing site and with helathy or COVID-infection phenotype
#seurat_subset <- seurat_object[,(seurat_object[[]]$batch=="Cambridge") & (seurat_object[[]]$Status=="Covid" | seurat_object[[]]$Status=="Healthy")]

# Select 20 random samples and subset to 500 random cells
#random_samples_ids <- sample(unique(seurat_subset[[]]$sample),20,replace=FALSE)
#seurat_subset <- seurat_subset[,seurat_subset[[]]$sample %in% random_samples_ids]
#cell_sample_tibble <- dplyr::tibble(sample_id = seurat_subset[[]]$sample, cell_barcode = rownames(seurat_subset[[]]))
#sample_subset <- cell_sample_tibble %>% dplyr::group_by(sample_id) %>% dplyr::slice_sample(n=500)
#seurat_subset <- seurat_subset[,rownames(seurat_subset[[]]) %in% sample_subset$cell_barcode]

# For unit-testing purposes we subsample to create a unit with 49 cells
#small_sample_id <- random_samples_ids[20]
#small_sample_drop_indices <- sample(which(seurat_subset[[]]$sample==small_sample_id),(500-49),replace=FALSE)
#small_sample_seurat <- seurat_subset[,-small_sample_drop_indices]

# we only keep the PCA embeddings of the 2000 selected HVGs
#example_object <- Seurat::DietSeurat(seurat_subset,counts=FALSE,data=TRUE,
#				     dimreducs="pca",misc=FALSE)
#testing_object <- Seurat::DietSeurat(small_sample_seurat,counts=FALSE,data=TRUE,
#					dimreducs="pca",misc=FALSE)

# The above function keeps the count data, which we remove manually
#example_object@assays <- list()
#testing_object@assays <- list()

# Finally, we save the object
#saveRDS(example_object,here::here("data","example_data.Rds"))
#saveRDS(testing_object,here::here("data","testing_data.Rds"))
