main_directory <- 'd:/master project/outcome/'
folders <- list.dirs(main_directory, recursive = FALSE)
final_result <- data.frame()
for (folder in folders) {
  if (basename(folder) == "ovary") {
    next
  }
setwd(folder)
names <- list.files(folder)
newantigens <- data.frame()
for (name in names) {
  epTable <- read.table(name,
                        sep='\t', stringsAsFactors = F, header=T,
                        col.names=c('Sample', 'LineID', 'Chrom', 'allelepos','REF',
                                    'ALTAll', 'FILTER', 'Gene', 'pos', 'hla', 'peptide', 'core', 'Of', 'Gp',
                                    'Gl', 'Ip', 'Il', 'Icore', 'Identity', 'Score','Affinity', 'Rank',
                                    'Cand','BindLevel', 'Novelty'),fill= T ,row.names=NULL)
  newantigens <- rbind(newantigens, subset(epTable, Novelty==1))
  newantigens.notdup <- newantigens[!duplicated(paste0(newantigens$Sample,newantigens$LineID)),]
}
if (nrow(newantigens) == 0) {
  next
}
sample_counts <- table(newantigens$Sample)
result_table <- data.frame(
  Sample = names(sample_counts),
  Total_Count = as.numeric(sample_counts),
  SB_Count = NA,
  SB_ratio = NA,
  WB_Count = NA,
  WB_ratio = NA,
  SB_nodup = NA
  
)
for (i in 1:nrow(result_table)) {
  sb_rows <- subset(newantigens, Sample == result_table$Sample[i] & BindLevel == "SB")
  wb_rows <- subset(newantigens, Sample == result_table$Sample[i] & BindLevel == "WB")
  sb_nodup <- subset(newantigens.notdup, Sample == result_table$Sample[i] & BindLevel == "SB")
  result_table$SB_Count[i] <- nrow(sb_rows)
  result_table$WB_Count[i] <- nrow(wb_rows)
  result_table$SB_nodup[i] <- nrow(sb_nodup)
}
result_table$SB_ratio <- result_table$SB_Count / result_table$Total_Count *100
result_table$WB_ratio <- result_table$WB_Count / result_table$Total_Count *100
result_table$Folder <- basename(folder)
final_result <- rbind(final_result, result_table)
}
library(ggplot2)
library(ggpubr)
final_result$Folder <- as.character(final_result$Folder)
add_significance <- function(plot, comparison, label, label_height) {
  comp_x <- sapply(comparison, function(folder) mean(which(final_result$Folder == folder)))
  comp_x <- mean(comp_x)
  plot + stat_compare_means(
    aes(label = ..p.signif..), 
    comparisons = list(comparison), 
    method = "t.test", 
    label.y = label_height,
    label.x = comp_x,
    inherit.aes = FALSE,
    label = label
  )
}
p <- ggplot(final_result, aes(x = Folder, y = SB_ratio)) +
  geom_violin(fill = "green", color = "blue") +
  geom_boxplot(width = 0.5, fill = "transparent", color = "blue") +
  labs(x = "Cancer", y = "SB Ratio(%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("SB Ratio")
p <- add_significance(p, c("brain", "breast"), "*", 60)
p <- add_significance(p, c("brain", "colon"), "**", 64)
p <- add_significance(p, c("brain", "lung"), "***", 68)
p <- add_significance(p, c("breast", "prostate"), "**", 72)
p <- add_significance(p, c("colon", "kidney"), "*", 60)
p <- add_significance(p, c("colon", "prostate"), "**", 64)
p <- add_significance(p, c("kidney", "lung"), "*", 60)
p <- add_significance(p, c("lung", "prostate"), "**", 60)
print(p)
main_directory <- 'd:/master project/VCF/'
vcf_folders <- list.dirs(main_directory, recursive = FALSE)
mutation_all <- data.frame(Sample = character(), Number = integer(), stringsAsFactors = FALSE)
for (vcf_folder in vcf_folders) {
  setwd(vcf_folder)
  vcf_names <- list.files(vcf_folder)
  for (vcf_name in vcf_names) {
    vcf_name_noext <- tools::file_path_sans_ext(vcf_name)
    vcf.df <- read.delim2(vcf_name, skip = 1, header = TRUE)
    vcf_num <- nrow(vcf.df)
    vcf_info <- data.frame(Sample = vcf_name_noext,Number=vcf_num,stringsAsFactors = FALSE)
    mutation_all <- rbind(mutation_all ,vcf_info)
  }
}
merged_data <- merge(final_result, mutation_all, by = "Sample")
merged_data$SB_Num_Ratio <- merged_data$SB_nodup/merged_data$Number *100
merged_data$Folder <- as.character(merged_data$Folder)
add_significance <- function(plot, comparison, label, label_height) {
  comp_x <- sapply(comparison, function(folder) mean(which(merged_data$Folder == folder)))
  comp_x <- mean(comp_x)
  plot + stat_compare_means(
    aes(label = ..p.signif..), 
    comparisons = list(comparison), 
    method = "t.test", 
    label.y = label_height,
    label.x = comp_x,
    inherit.aes = FALSE,
    label = label
  )
}
g <- ggplot(merged_data, aes(x = Folder, y = SB_Num_Ratio)) +
  geom_violin(fill = "green", color = "blue") +
  geom_boxplot(width = 0.1, fill = "transparent", color = "blue") +
  labs(x = "Cancer", y = "SB/Mutation Ratio(%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ggtitle("SB/Mutation Ratio")

g <- add_significance(g, c("breast", "colon"), "*", 60)
g <- add_significance(g, c("kidney", "lung"), "*", 60)
g <- add_significance(g, c("colon", "lung"), "**", 65)
print(g)
