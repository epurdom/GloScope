data("example_data")

set.seed(1)
library(dplyr)
library(stringr)
sub_pat <- sample(unique(example_data$donor_label),5)

sub_data <- example_data[example_data$donor_label %in% sub_pat,]
sub_meta <- unique(sub_data[, c("donor_label", "joint_region_label")])

sample_name = as.character(unique(sub_meta[, "donor_label"]))
sub_meta[,"donor_label"] = as.character(sub_meta[,"donor_label"])
df_list = split(sub_data, sub_data[,"donor_label"])
df_list = lapply(df_list, function(y) y[,str_detect(colnames(y), "PC")])
df_list = lapply(df_list, function(y) as.matrix(y[,1:10]))
