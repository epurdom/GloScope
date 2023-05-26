data(example_data)
set.seed(2)

sample_ids <- example_data$metadata$sample_id
subsample_patients <- sample(unique(sample_ids),3)
subsample_metadata <- example_data$metadata[sample_ids %in% subsample_patients,]
subsample_data <- example_data$pca_embeddings[sample_ids %in% subsample_patients,]
subsample_data_subset <- subsample_data[,1:10] # Pick the first 10 PCs

undersized_patient <- subsample_patients[1]
sample_drop_indices <- sample(which(subsample_metadata==undersized_patient),(500-49),replace=FALSE)
undersized_metadata <- subsample_metadata[-sample_drop_indices,]
undersized_data <- subsample_data[-sample_drop_indices,]
undersized_data_subset <- undersized_data[,1:10]
