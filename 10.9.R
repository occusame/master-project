read_maf_from_gz <- function(gz_file) {
  sam <- gzfile(gz_file, "rt")
  maf_text <- read.delim(sam,comment.char = '#')
  close(sam)
  return(maf_text)
}

setwd('d:/master project/portal samples data/brain/')
current_directory <- getwd()
names <- list.files(current_directory, full.names = TRUE)
names <- names[!grepl('antigen_presenting_genes|outcome', names)]
num <- length(names)-1

for (i in 1:num) {
  folder <- names[i]
  gz_files <- list.files(folder, pattern = "\\.gz$", full.names = TRUE)
  for (gz_file in gz_files) {
    maf_text <- read_maf_from_gz(gz_file)
    df_name <- paste("case",i,".df", sep = "_")
    assign(df_name, maf_text)
    cat("Processing is number" , i , "floder", gz_file, "\n")
  }
}

gene.df <- read.delim('antigen_presenting_genes.txt')
genes <- gene.df$HGNC.symbol
mut_escape_list <- list()
mut_numbers <- matrix(NA, nrow = num, ncol = 3)
names(mut_numbers) <- c('Patient','EscapeNumber','EscapeGenes')
gene.df <- read.delim('antigen_presenting_genes.txt')
genes <- gene.df$HGNC.symbol
mut_escape_list <- list()

for (i in 1:num) {
  df_name <- paste("case", i, ".df", sep = "_")
  mut_file.df <- get(df_name)
  table(mut_file.df$Variant_Classification)
  table(mut_file.df$IMPACT)
  mut_file.escape.df <- subset(mut_file.df, (Hugo_Symbol %in% genes) & !(Variant_Classification %in% c('Silent','Intron')))
  mut_escape_list[[i]] <- mut_file.escape.df[,c(1,4,9,16)]
  escape_number <- nrow(mut_file.escape.df)
  escape_genes <- paste0(mut_file.escape.df$Hugo_Symbol, collapse=';')
  mut_numbers[i,] <- c(df_name,escape_number, escape_genes)
}

for (i in 1:num) {
  cat("Results for iteration", i, ":\n")
  print(mut_escape_list[[i]])
}

existing_file <- 'outcome.csv'
mut_escape_short_list <- list()
for (i in 1:num){
  x <- mut_escape_list[[i]]
  mut_escape_short_list[[i]] <- x[,c('Hugo_Symbol','Tumor_Sample_Barcode')]
}

combined_mutations <- do.call(rbind, mut_escape_short_list)
new_data_frame <- mut_escape_list
write.csv(combined_mutations, file = existing_file, row.names = FALSE)
