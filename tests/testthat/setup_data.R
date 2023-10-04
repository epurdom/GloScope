data(example_SCE)

sample_ids <- SingleCellExperiment::colData(example_SCE)$sample_id
set.seed(2)
subsample_patients <- sample(unique(sample_ids),4)
subsample_metadata <- as.data.frame(SingleCellExperiment::colData(example_SCE)[sample_ids %in% subsample_patients,])
subsample_data <- SingleCellExperiment::reducedDim(example_SCE,"PCA")[sample_ids %in% subsample_patients,]
subsample_data_subset <- subsample_data[,1:10] # Pick the first 10 PCs

undersized_patient <- subsample_patients[1]
set.seed(2)
sample_drop_indices <- sample(which(subsample_metadata==undersized_patient),(500-49),replace=FALSE)
undersized_metadata <- subsample_metadata[-sample_drop_indices,]
undersized_data <- subsample_data[-sample_drop_indices,]
undersized_data_subset <- undersized_data[,1:10]
