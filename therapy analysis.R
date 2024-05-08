gene.df <- read.delim('d:/master project/antigen_presenting_genes.txt')
genes <- gene.df$HGNC.symbol
main_directory <- 'd:/master project/cbioportal/'
setwd(main_directory)
fold_names <- list.files(main_directory)
fold_names <- fold_names[!grepl('.tar.gz|lung|Colorectal Adenocarcinoma|therapy information', fold_names)]
result_data <- data.frame()

for (fold_name in fold_names) {
  mutation_file <- file.path(fold_name, 'data_mutations.txt')
  mut.file <- read.delim(mutation_file, comment.char = '#')
  
  gene_column <- which(colnames(mut.file) %in% "Hugo_Symbol")
  variant_column <- which(colnames(mut.file) %in% "Variant_Classification")
  barcode_column <- which(colnames(mut.file) %in% "Tumor_Sample_Barcode")
  
  mut.escape.df <- subset(mut.file, (Hugo_Symbol %in% genes)& !(Variant_Classification %in% c('Silent','Intron')))
  mut_escape <- mut.escape.df[,c(barcode_column,gene_column, variant_column)]
  escape_number <- nrow(mut.escape.df)
  escape_genes <- paste0(mut.escape.df$Hugo_Symbol, collapse=';')
  print(mut_escape)
  allcount <- nrow(mut.file)
  print(allcount)
  mutcount <- table(mut_escape$Hugo_Symbol)
  count <- numeric(length = length(genes)); names(count) = genes
  if (length(mutcount) > 0) {
    count[match(names(mutcount), names(count))] <- mutcount / allcount * 100
  }
  data_to_plot <- data.frame(
    Hugo_Symbol = names(count),
    Percentage = as.numeric(count),
    Source = fold_name
  )
  result_data <- rbind(result_data, data_to_plot)
}
print(result_data)