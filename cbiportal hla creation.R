hla.df <- read.delim('d:/master project/cbiportal vcf/cbiportal_HLA_data.txt', header = FALSE)
# 计算需要创建的文件数量
file_count <- ceiling(nrow(hla.df) / 30)

# 循环遍历每个文件
for (i in 1:file_count) {
  # 计算当前子集的起始和结束行号
  start_row <- (i - 1) * 30 + 1
  end_row <- min(i * 30, nrow(hla.df))
  
  # 提取子集
  subset_df <- hla.df[start_row:end_row, ]
  
  # 构建文件名
  file_name <- sprintf("d:/master project/cbiportal hla/hla_data_%02d.txt", i)
  
  # 保存文件，不包括列名称
  write.table(subset_df, file_name, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
}
