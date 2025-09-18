data(example_SCE)

sample_ids <- SingleCellExperiment::colData(example_SCE)$sample_id
set.seed(2)
subsample_patients <- sample(unique(sample_ids),4)
subsample_metadata <- as.data.frame(SingleCellExperiment::colData(example_SCE)[sample_ids %in% subsample_patients,])
subsample_data <- SingleCellExperiment::reducedDim(example_SCE,"PCA")[sample_ids %in% subsample_patients,]
subsample_data_subset <- subsample_data[,1:10] # Pick the first 10 PCs

undersized_patient <- subsample_patients[1]
set.seed(2)
sample_drop_indices <- sample(which(subsample_metadata$sample_id==undersized_patient),(500-49),replace=FALSE)
undersized_metadata <- subsample_metadata[-sample_drop_indices,]
undersized_data <- subsample_data[-sample_drop_indices,]
undersized_data_subset <- undersized_data[,1:10]

#for patient level tests need a gloscope matrix to test on
expect_silent(dist_mat <- gloscope(embedding_matrix=subsample_data_subset, cell_sample_ids=subsample_metadata$sample_id,dens="KNN"))
pat_info <- unique(subsample_metadata[,c(1,2)])
pat_info$group <- c("A","A","B","B")

pat_info_full<-as.data.frame(unique(SingleCellExperiment::colData(example_SCE)[,c("sample_id","phenotype")]))
pat_info_full$group<-rep(c("A","B"),each=10)
#create fake dist_mat of 20x20
n<-nrow(pat_info_full)
rvalues<-round(abs(rnorm(n*(n-1)/2)),3)
dist_mat_full<-matrix(0,nrow=n,ncol=n)
dist_mat_full[upper.tri(dist_mat_full)]<-rvalues
dist_mat_full[lower.tri(dist_mat_full)]<-t(dist_mat_full)[lower.tri(dist_mat_full)]
rownames(dist_mat_full)<-colnames(dist_mat_full)<-pat_info_full$sample_id