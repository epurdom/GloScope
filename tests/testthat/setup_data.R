data(example_data)
set.seed(2)

sub_pat <- sample(unique(example_data$patient_id),3)

sub_data <- example_data[example_data$patient_id %in% sub_pat,]
#sub_meta <- unique(sub_data[, c("patient_id", "Status")])

#sample_name = as.character(unique(sub_meta[, "patient_id"]))
#sub_meta[,"patient_id"] = as.character(sub_meta[,"patient_id"])
#df_list = split(sub_data, sub_data[,"patient_id"])
#df_list = lapply(df_list, function(y) y[,str_detect(colnames(y), "PC")])
#df_list = lapply(df_list, function(y) as.matrix(y[,1:10]))
##################################################

sub_ids_40 <- sample(which(sub_data$patient_id == sub_pat[1]),40)
sub_data_40 <- rbind(sub_data[sub_data$patient_id!=sub_pat[1],],
                       sub_data[sub_ids_40,])

sub_ids_350 <- sample(which(sub_data$patient_id == sub_pat[1]),350)
sub_data_350 <- rbind(sub_data[sub_data$patient_id!=sub_pat[1],],
                     sub_data[sub_ids_350,])
