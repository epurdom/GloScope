# The examples, test cases, and vignette for this package use a subset of
# peripheral blood mononuclear cells (PBMCs) scRNA-Seq data from
# COVID-infected patients and healthy controls presented in
# "Single-cell multi-omics analysis of the immune response in COVID-19"
# by Stephenson et al. (2021) [doi: 10.1038/s41591-021-01329-2]. This 
# script recreates the two .Rda files found in the `data/` folder.

#\dontrun{
# An `.h5ad` file containing the processed counts in Seurat format is
# downloaded from the European Bioinformatics Institute with
# Accession E-MTAB-10026. Note that this file is on the order of 8 GB
# and can take several hours to download via FTP.

here::i_am("inst/extdata/create_example_data.R")
library(SingleCellExperiment)
library(zellkonverter)

# Execute a shell script provided by the European Bioinformatics Institute
system("bash inst/extdata/ftp_data.sh")

# Load .h5ad file as SingleCellExperiment
full_sce <- zellkonverter::readH5AD(here::here("inst","extdata","covid_portal_210320_with_raw.h5ad"))

# Subset to 500 cells from 20 patients
subset_cell_labels <- readRDS(here::here("inst","extdata","subset_cell_labels.Rds"))
subset_sce <- full_sce[,subset_cell_labels]

# Create a minimal SCE with only metadata and PCA embeddings
subset_metadata <- colData(subset_sce)[,c("patient_id","Status")]
colnames(subset_metadata) <- c("sample_id","phenotype")
subset_metadata[] <- lapply(subset_metadata, as.character)

subset_pca <- SingleCellExperiment::reducedDim(subset_sce,"X_pca")
pca_colnames <- paste0(rep("PC_",50),1:50)
colnames(subset_pca) <- pca_colnames
zero_counts <- matrix(0,2,10000) # Necessary to add PCA to SCE
colnames(zero_counts) <- subset_cell_labels

example_SCE <- SingleCellExperiment::SingleCellExperiment(assays=list(zero_counts), reducedDims=list(PCA=subset_pca))
colData(example_SCE) <- subset_metadata
example_SCE@assays@data@listData <- list() # Remove zero count matrix

save(example_SCE,file=here::here("inst","extdata","example_SCE.rda"))
#}
