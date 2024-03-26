gene.df <- read.delim('d:/master project/antigen_presenting_genes.txt')
genes <- gene.df$HGNC.symbol
main_directory <- 'd:/master project/cbioportal/lung/'
setwd(main_directory)
main_directory <- getwd()
fold_names <- list.files(main_directory, full.names = FALSE)
fold_names <- fold_names[!grepl('.tar.gz', fold_names)]
combined_data <- data.frame()
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
  percent <- mutcount/allcount *100
  print(percent)
  plot_data <- data.frame(Gene = names(percent), Percentage = as.numeric(percent), MutCount = as.numeric(mutcount),Source = fold_name)
  combined_data <- rbind(combined_data, plot_data)
}

library(ggplot2)
ggplot(combined_data, aes(x = Gene, y = Percentage, fill = factor(Source))) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = MutCount), 
            position = position_dodge(width = 0.9), 
            vjust = -0.5, size = 2) +
  labs(x = "Hugo_Symbol", y = "Percentage") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5))

