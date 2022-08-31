data("example_data")

set.seed(1)
library(dplyr)
library(stringr)
sub_pat <- sample(unique(example_data$patient_id),5)

sub_data <- example_data[example_data$patient_id %in% sub_pat,]
sub_meta <- unique(sub_data[, c("patient_id", "Status")])

sample_name = as.character(unique(sub_meta[, "patient_id"]))
sub_meta[,"patient_id"] = as.character(sub_meta[,"patient_id"])
df_list = split(sub_data, sub_data[,"patient_id"])
df_list = lapply(df_list, function(y) y[,str_detect(colnames(y), "PC")])
df_list = lapply(df_list, function(y) as.matrix(y[,1:10]))
