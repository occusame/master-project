hla.df <- read.delim('d:/master project/cbiportal vcf/cbiportal_HLA_data.txt', header = FALSE)

file_count <- ceiling(nrow(hla.df) / 30)

for (i in 1:file_count) {
  start_row <- (i - 1) * 30 + 1
  end_row <- min(i * 30, nrow(hla.df))
  
  subset_df <- hla.df[start_row:end_row, ]
  
  file_name <- sprintf("d:/master project/cbiportal hla/hla_data_%02d.txt", i)
  
  write.table(subset_df, file_name, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
}
