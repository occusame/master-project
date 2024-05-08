main_directory <- 'd:/master project/cbioportal'
fold_names <- c()
fold_names <- c(fold_names, list.files(main_directory, full.names = TRUE))
fold_names <- fold_names[!grepl('lung|Colorectal Adenocarcinoma', fold_names)]
fold_names <- c(fold_names, list.files(file.path(main_directory, "lung"), full.names = TRUE))
fold_names <- c(fold_names, list.files(file.path(main_directory, "Colorectal Adenocarcinoma"), full.names = TRUE))
fold_names <- fold_names[!grepl('.tar.gz|.txt', fold_names)]

names <- c()
for (fold_name in fold_names) {
  mutation_file <- file.path(fold_name, 'data_mutations.txt')
  mut.file <- read.delim(mutation_file, comment.char = '#')
  if ("t_depth" %in% names(mut.file)) {
    if ("dbSNP_RS" %in% names(mut.file)) {
      select <- mut.file[,c('Tumor_Sample_Barcode','Chromosome','Start_Position','dbSNP_RS','Reference_Allele','Tumor_Seq_Allele2','t_depth','t_ref_count','t_alt_count')]
    } else {
      select <- mut.file[,c('Tumor_Sample_Barcode','Chromosome','Start_Position','Reference_Allele','Tumor_Seq_Allele2','t_depth','t_ref_count','t_alt_count')]
    }
    select$INFO <- paste0(select$t_depth, ':', select$t_ref_count, ':', select$t_alt_count)
    select <- select[, !colnames(select) %in% c('t_depth', 't_ref_count', 't_alt_count')]
  } else {
    select <- mut.file[,c('Tumor_Sample_Barcode','Chromosome','Start_Position','dbSNP_RS','Reference_Allele','Tumor_Seq_Allele2','t_ref_count','t_alt_count')]
    select$INFO <- paste0(select$t_ref_count,':',select$t_alt_count)
    select <- select[,!colnames(select) %in% c('t_ref_count','t_alt_count')]
  }
  select <- cbind(select[, 1:6], QUAL = '50', FILTER = 'PASS', INFO = select$INFO, FORMAT = 'GT:', barcode = '0/1:')
  select$FORMAT <- paste0(select$FORMAT, select$Tumor_Seq_Allele2)
  select$barcode <- paste0(select$barcode, mut.file$Tumor_Seq_Allele2)
  splits <- split(select, select$Tumor_Sample_Barcode)
  for (barcode in names(splits)) {
    split_data <- splits[[barcode]]
    names(split_data)[names(split_data) == "barcode"] <- barcode
    split_data <- split_data[,-1]
    names(split_data)[1:5] <- c('#CHROM', 'POS', 'ID', 'REF', 'ALT')
    output_folder <- 'd:/master project/cbiportal vcf/'
    output_file <- file.path(output_folder, paste0(barcode, ".vcf"))
    writeLines('##fileformat=VCFv4.2', output_file)
    write.table(split_data, file = output_file, sep = '\t', quote = FALSE, col.names = TRUE, row.names = FALSE, append = TRUE)
    hla_data <- data.frame(Tumor_Sample_Barcode = barcode,
                           HLA_A_02_01 = rep('HLA_a_02_01', length(names)),
                           HLA_A_01_01 = rep('HLA_A_01_01', length(names)),
                           HLA_C_01_02 = rep('HLA_C_01_02', length(names)),
                           HLA_C_07_01 = rep('HLA_C_07_01', length(names)),
                           HLA_B_08_01 = rep('HLA_B_08_01', length(names)),
                           HLA_B_07_02 = rep('HLA_B_07_02', length(names)),
                           stringsAsFactors = FALSE)
    hla_output_folder <- 'd:/master project/cbiportal hla detailed/'
    cbiportal <- file.path(hla_output_folder, "cbiportal_",barcode,"_HLA_data.txt")
    write.table(hla_data, file = cbiportal, sep = '\t', quote = FALSE, col.names = FALSE, row.names = FALSE)
  }
}


folder_path <- 'd:/master project/cbiportal hla detailed/'
txt_files <- list.files(folder_path, pattern = "cbiportal", full.names = TRUE)
files_merge <- 20
merge_number <- ceiling(length(txt_files) / files_merge)
for (i in 1:merge_number) {
  start_index <- (i - 1) * files_merge + 1
  end_index <- min(i * files_merge, length(txt_files))
  files_to_merge <- txt_files[start_index:end_index]
  combined_data <- lapply(files_to_merge, readLines)
  combined_data <- unlist(combined_data)
  output_file <- file.path('d:/master project/cbiportal hla/', paste0("cbiportal_", i, ".txt"))
  writeLines(combined_data, con = output_file)
  cat("Group", i, "merged successfully.\n")
}