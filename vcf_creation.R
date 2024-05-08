sample.maf <- read.delim('d:/master project/maf sample.maf', comment.char = '#')
select <- sample.maf[,c('Chromosome','Start_Position','dbSNP_RS','Reference_Allele','Allele','t_depth','t_ref_count','t_alt_count')]
select$INFO <- paste0(select$t_depth,':',select$t_ref_count,':',select$t_alt_count)
select <- select[,!colnames(select) %in% c('t_depth','t_ref_count','t_alt_count')]
select <- cbind(select[, 1:5], QUAL = '50', FILTER = 'PASS', INFO = select$INFO, FORMAT = 'GT:', `TCGA-S9-A7R4-01A-12D-A34J-08` = '0/1:')
select$FORMAT <- paste0(select$FORMAT, select$Allele)
select$`TCGA-S9-A7R4-01A-12D-A34J-08` <- paste0(select$`TCGA-S9-A7R4-01A-12D-A34J-08`, select$Allele)
names(select)[1:5] <- c('#CHROM','POS','ID','REF','ALT')
writeLines('##fileformat=VCFv4.2', 'd:/master project/sample.vcf')
write.table(select,file = 'd:/master project/sample.vcf',sep = '\t',quote = FALSE,col.names = TRUE,row.names = FALSE,append = TRUE)

sa.df <-read.delim('d:/master project/sample.vcf', comment.char = '#')
outcome.df <- read.delim2('d:/master project/sample.neoantigens.txt')
