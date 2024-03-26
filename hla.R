list.df <- read.csv('d:/master project/OptiTypeCallsHLA_20171207.tsv')
library(dplyr)
library(stringr)
list.df <- list.df %>%
  mutate_all(~str_replace_all(., "A\\*", "HLA_a_")) %>%
  mutate_all(~str_replace_all(., "B\\*", "HLA_b_")) %>%
  mutate_all(~str_replace_all(., "C\\*", "HLA_c_")) %>%
  mutate_all(~str_replace_all(., ":", "_")) %>%
  select(-Reads, -Objective)
list.df <- list.df[, c("aliquot_id", setdiff(names(list.df), "aliquot_id"))]
list.df$aliquot_id <- substr(list.df$aliquot_id, 1, 12)
read_maf_from_gz <- function(gz_file) {
  sam <- gzfile(gz_file, "rt")
  maf_text <- read.delim(sam,comment.char = '#')
  close(sam)
  return(maf_text)
}
main_directory <- 'd:/master project/portal samples data/prostate gland/'
setwd(main_directory)
hla.df <- data.frame()
names <- list.files(main_directory, full.names = TRUE)
names <- names[!grepl('.txt|outcome', names)]
num <- length(names)
for (i in 1:num) {
  folder <- names[i]
  gz_files <- list.files(folder, pattern = "\\.gz$", full.names = TRUE)
  for (gz_file in gz_files) {
    maf_text <- read_maf_from_gz(gz_file)
    maf_text$Tumor_Sample_Barcode <- substr(maf_text$Tumor_Sample_Barcode, 1, 12)
    df_name <- paste("case",i,".df", sep = "_")
    assign(df_name, maf_text)
    select <- get(df_name)[, c('Chromosome', 'Start_Position', 'dbSNP_RS', 'Reference_Allele', 'Allele', 't_depth', 't_ref_count', 't_alt_count')]
    if (nrow(select) == 0) {
      warning(paste("Empty data frame for case", i))
      next
    }
    select$INFO <- paste0(select$t_depth, ':', select$t_ref_count, ':', select$t_alt_count)
    select <- select[, !colnames(select) %in% c('t_depth', 't_ref_count', 't_alt_count')]
    select <- cbind(select[, 1:5], QUAL = '50', FILTER = 'PASS', INFO = select$INFO, FORMAT = 'GT:', barcode = '0/1:')
    hla.df <- list.df[list.df$aliquot_id %in% maf_text$Tumor_Sample_Barcode,]
    hla.df <- hla.df[!duplicated(hla.df$aliquot_id),]
    hlaoutput <- file.path('d:/master project/HLA detailed/', paste("prostate_hla_case",i,".txt", sep = "_"))
    write.table(hla.df,file = hlaoutput,sep = '\t',quote = FALSE,col.names = FALSE,row.names = FALSE,append = TRUE)
    select$FORMAT <- paste0(select$FORMAT, select$Allele)
    select$barcode <- paste0(select$barcode, select$Allele)
    barcode_name <- as.character(get(df_name)$Tumor_Sample_Barcode)
    barcode_name <- barcode_name[!duplicated(barcode_name)]
    barcode_name <- gsub("/$", "", barcode_name)
    names(select)[which(names(select) == "barcode")] <- barcode_name
    names(select)[1:5] <- c('#CHROM', 'POS', 'ID', 'REF', 'ALT')
    output_folder <- 'd:/master project/VCF/prostate'
    output_file <- file.path(output_folder, paste0(barcode_name, ".vcf"))
    writeLines('##fileformat=VCFv4.2', output_file)
    write.table(select, file = output_file, sep = '\t', quote = FALSE, col.names = TRUE, row.names = FALSE, append = TRUE)
  }
}

folder_path <- 'd:/master project/HLA detailed/'
txt_files <- list.files(folder_path, pattern = "prostate_hla_case", full.names = TRUE)
files_merge <- 20
merge_number <- ceiling(length(txt_files) / files_merge)
for (i in 1:merge_number) {
  start_index <- (i - 1) * files_merge + 1
  end_index <- min(i * files_merge, length(txt_files))
  files_to_merge <- txt_files[start_index:end_index]
  combined_data <- lapply(files_to_merge, readLines)
  combined_data <- unlist(combined_data)
  output_file <- file.path('d:/master project/HLA/prostate', paste0("HLA_prostate_", i, ".txt"))
  writeLines(combined_data, con = output_file)
  cat("Group", i, "merged successfully.\n")
}
