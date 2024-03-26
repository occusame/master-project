gene.df <- read.delim('d:/master project/antigen_presenting_genes.txt')
genes <- gene.df$HGNC.symbol
setwd('d:/master project/cbioportal/Colorectal Adenocarcinoma/coadread_dfci_2016/')
samples.df <- read.delim('data_clinical_sample.txt')
samples.df <- samples.df[-c(1:4), ]
patient.df <- read.delim('data_clinical_patient.txt')
patient.df <- patient.df[-c(1:4), ]

mut.df <- read.delim('data_mutations.txt', comment.char = '#')
mut.escape.df <- subset(mut.df, (Hugo_Symbol %in% genes)& !(Variant_Classification %in% c('Silent','Intron')))
gene.column <- which(colnames(mut.df) %in% "Hugo_Symbol")
variant.column <- which(colnames(mut.df) %in% "Variant_Classification")
barcode.column <- which(colnames(mut.df) %in% "Tumor_Sample_Barcode")
mut.escape <- mut.escape.df[,c(barcode.column,gene.column, variant.column)]

infor.columns <- samples.df[,c('X.Sample.Identifier','Oncotree.Code','Cancer.Type')]
add.col <- patient.df[,c('X.Patient.Id','Adjuvant.Treatment')]

total.df <- merge(infor.columns, add.col, by.x = "X.Patient.Identifier", by.y = "X.Patient.Id", all.x = TRUE)
total.df <- total.df[, -which(names(total.df) == "X.Patient.Identifier")]
colnames(infor.columns)[colnames(infor.columns) == "Sample.Type"] <- "Therapy"
infor.columns$Therapy[infor.columns$Therapy == 'Primary'] <- 'NTX'
total.df$Therapy[!(total.df$Therapy %in% c('NTX'))] <- 'TX'

part.df<- merge( mut.escape, infor.columns, by.x = "Tumor_Sample_Barcode", by.y = "X.Sample.Identifier", all.x = TRUE)
write.csv(part.df, file = "information.csv", row.names = FALSE)
write.csv(infor.columns, file = "infor.csv", row.names = FALSE)
