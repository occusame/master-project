setwd('d:/master project/cbioportal/lung/sclc_ucologne_2015/')
samples.df <- read.delim('data_clinical_sample.txt')
samples.df <- samples.df[-c(1:4), ]
patient.df <- read.delim('data_clinical_patient.txt')
patient.df <- patient.df[-c(1:4), ]

mut.df <- read.delim('data_mutations.txt', comment.char = '#')
gene.df <- read.delim('d:/master project/antigen_presenting_genes.txt')
genes <- gene.df$HGNC.symbol
mut.escape.df <- subset(mut.df, (Hugo_Symbol %in% genes)& !(Variant_Classification %in% c('Silent','Intron')))
gene.column <- which(colnames(mut.df) %in% "Hugo_Symbol")
variant.column <- which(colnames(mut.df) %in% "Variant_Classification")
barcode.column <- which(colnames(mut.df) %in% "Tumor_Sample_Barcode")
mut.escape <- mut.escape.df[,c(barcode.column,gene.column, variant.column)]

infor.columns <- samples.df[,c('Sample.Identifier','Oncotree.Code','Cancer.Type')]
add.col <- patient.df[,c('X.Patient.Identifier','Previous.Treatment')]
colnames(merged_data)[colnames(merged_data) == "Previous.Treatment"] <- "Therapy"
infor.columns$Therapy <- 'NTX'
infor.columns$Therapy[infor.columns$Therapy == 'Primary'] <- 'NTX'
merged_data$Therapy[!(merged_data$Therapy %in% c('TX', 'NTX'))] <- 'Unknown'



merged_data <- merge(mut.escape, infor.columns, by.x = "Tumor_Sample_Barcode", by.y = "Sample.Identifier", all.x = TRUE)
merged_data <- merge(merged_data, infor2.columns, by.x = "Tumor_Sample_Barcode", by.y = "X.Patient.Identifier", all.x = TRUE)
write.csv(merged_data, file = "information.csv", row.names = FALSE)
write.csv(infor.columns, file = "infor.csv", row.names = FALSE)
