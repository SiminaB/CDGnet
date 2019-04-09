library(tidyverse)
library(readr)

products <- read_delim("data/Products.txt", 
                       "\t", escape_double = FALSE,
                       col_types = cols(ApplNo = col_integer(), ProductNo = col_integer()), 
                       trim_ws = TRUE)

submissions <- read_delim("data/Submissions.txt", 
                          "\t", escape_double = FALSE, 
                          col_types = cols(ApplNo = col_integer(), SubmissionNo = col_integer()), 
                          trim_ws = TRUE)


merge_data <- merge(products, submissions, by="ApplNo")
merge_data <- as.data.frame(merge_data)

merge_data$drugId <- paste0(merge_data$ApplNo, "-", merge_data$ProductNo)

all_approved_drugs <- filter(merge_data, merge_data$SubmissionStatus %in% c("AP", "TA") )

# all approved drugs from FDA
drug_names <- unique(all_approved_drugs$DrugName)


# some of these drugs might be discontinued. filtering those out

marketing_status <- read_delim("data/MarketingStatus.txt", 
                               "\t", escape_double = FALSE, 
                               col_types = cols(ApplNo = col_integer(), ProductNo = col_integer()), 
                               trim_ws = TRUE)

marketing_status$drugId <- paste0(marketing_status$ApplNo, "-", marketing_status$ProductNo)


filter_discontinued <- filter(marketing_status, marketing_status$MarketingStatusID == 3)

all_approved_drugs_filter <- filter(merge_data, !merge_data$drugId %in% filter_discontinued$drugId)

drug_names_filtered <- unique(all_approved_drugs_filter$DrugName)

## Orange Book list

olist_drugs <- read_delim("orange_book_list/products.txt", 
                       "~", escape_double = FALSE, trim_ws = TRUE)

olist_drug <- unique(olist_drugs$Trade_Name)

# write to file

write.table(as.data.frame(drug_names),file="drug_names.txt", quote=F,sep=",",row.names=F)

write.table(as.data.frame(drug_names_filtered),file="drug_names_filtered.txt", quote=F,sep=",",row.names=F)

write.table(as.data.frame(olist_drug),file="olist_drug_names.txt", quote=F,sep=",",row.names=F)
