treatment <- data.frame()
setwd('d:/master project/cbioportal/')
blca_tcga_pub_2017.df <- read.delim('blca_tcga_pub_2017/data_mutations.txt', comment.char = '#')
blca_2017 <- data.frame(sample = blca_tcga_pub_2017.df$Tumor_Sample_Barcode)
blca_2017$treatment <- "NTX"
blca_2017$project <- 'blca_tcga_pub_2017'
treatment <- rbind(treatment,blca_2017)
brca_pareja_msk_2020.df <- read.delim('brca_pareja_msk_2020/data_mutations.txt', comment.char = '#')
brca_2020 <- data.frame(sample = brca_pareja_msk_2020.df$Tumor_Sample_Barcode)
brca_2020$treatment <- "NTX"
brca_2020$project <- 'brca_pareja_msk_2020'
treatment <- rbind(treatment,brca_2020)
ccrcc_dfci_2019.df <- read.delim('ccrcc_dfci_2019/data_mutations.txt', comment.char = '#')
ccrcc_2019 <- data.frame(sample = ccrcc_dfci_2019.df$Tumor_Sample_Barcode)
ccrcc_2019$treatment <- "TX"
ccrcc_2019$project <- 'ccrcc_dfci_2019'
treatment <- rbind(treatment,ccrcc_2019)
cll_broad_2015.df <- read.delim('cll_broad_2015/data_mutations.txt', comment.char = '#')
cllsample <- read.delim('cll_broad_2015/data_clinical_sample.txt', comment.char = '#')
cllpatient <- read.delim('cll_broad_2015/data_clinical_patient.txt', comment.char = '#')
cll_2015 <-  unique(data.frame(sample = cll_broad_2015.df$Tumor_Sample_Barcode))
merged_df <- merge(cll_2015, cllsample, by.x = "sample", by.y = "SAMPLE_ID", all.x = TRUE)
merged_df <- merge(merged_df,cllpatient, by.x = "PATIENT_ID", by.y = "PATIENT_ID", all.x = TRUE)
cll_2015 <- merged_df[,c('sample','PRIOR_TREATMENT')]
colnames(cll_2015)[colnames(cll_2015) == "PRIOR_TREATMENT"] <- "treatment"
cll_2015$project <- 'cll_broad_2015'
cll_2015$treatment[cll_2015$treatment == 'prior therapy'] <- "TX"
cll_2015$treatment[cll_2015$treatment == 'treatment naive'] <- "TX"
treatment <- rbind(treatment,cll_2015)
gbc_shanghai_2014.df <- read.delim('gbc_shanghai_2014/data_mutations.txt', comment.char = '#')
gbc_2014 <-  data.frame(sample = gbc_shanghai_2014.df$Tumor_Sample_Barcode)
gbc_2014$treatment <- "TX"
gbc_2014$project <- 'gbc_shanghai_2014'
treatment <- rbind(treatment,gbc_2014)
gbm_cptac_2021.df <- read.delim('gbm_cptac_2021/data_mutations.txt', comment.char = '#')
gbm <- read.delim('gbm_cptac_2021/data_clinical_patient.txt', comment.char = '#')
gbm_cptac_2021<-  data.frame(sample = gbm_cptac_2021.df$Tumor_Sample_Barcode)
gbm_subset <- gbm[, c("PATIENT_ID", "TUMOR_REOCCUR_AFTER_TREATMENT")]
gbm_cptac_2021<- merge(gbm_cptac_2021, gbm_subset, by.x = "sample", by.y = "PATIENT_ID", all.x = TRUE)
names(gbm_cptac_2021)[names(gbm_cptac_2021) == "TUMOR_REOCCUR_AFTER_TREATMENT"] <- "treatment"
gbm_cptac_2021$treatment[gbm_cptac_2021$treatment == FALSE] <- "TX"
gbm_cptac_2021$treatment[gbm_cptac_2021$treatment == TRUE] <- "TX"
gbm_cptac_2021$treatment[is.na(gbm_cptac_2021$treatment)] <- "NTX" 
gbm_cptac_2021$project <- 'gbm_cptac_2021'
treatment <- rbind(treatment,gbm_cptac_2021)
hcc_meric_2021.df <- read.delim('hcc_meric_2021/data_mutations.txt', comment.char = '#')
hcc_meric_2021 <- data.frame(sample = hcc_meric_2021.df$Tumor_Sample_Barcode)
hcc_meric_2021$treatment <- "NTX"
hcc_meric_2021$project <- 'hcc_meric_2021'
treatment <- rbind(treatment,hcc_meric_2021)
pptc_2019.df <- read.delim('pptc_2019/data_mutations.txt', comment.char = '#')
pptc_2019 <- data.frame(sample = pptc_2019.df$Tumor_Sample_Barcode)
unique_rows <- !duplicated(treatment$sample)
treatment_unique <- treatment[unique_rows, ]
write.csv(treatment_unique, file = "therapy_information.csv", row.names = FALSE)
