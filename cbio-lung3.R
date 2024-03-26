main_directory <- 'd:/master project/cbioportal/lung/'
setwd(main_directory)
main_directory <- getwd()
fold_names <- list.files(main_directory, full.names = FALSE)
fold_names <- fold_names[!grepl('.tar.gz|txt', fold_names)]
part.df <- data.frame()
total.df <- data.frame()
information1.df <- read.csv('d:/master project/cbioportal/lung/luad_broad/information.csv')
infor1.df <- read.csv('d:/master project/cbioportal/lung/luad_broad/infor.csv')
for (fold_name in fold_names) {
  information_file <- file.path(fold_name, "information.csv")
  information_data <- read.csv(information_file)
  infor_file <- file.path(fold_name, "infor.csv")
  infor_data <- read.csv(infor_file)
  colnames(information_data) <- colnames(information1.df)
  colnames(infor_data) <- colnames(infor1.df)
  part.df <- rbind(part.df,information_data)
  total.df <- rbind(total.df, infor_data)
}
part.df$Therapy[part.df$Therapy == 'No Therapy'] <- 'NTX'
part.df$Therapy[part.df$Therapy == 'Untreated'] <- 'NTX'
part.df$Therapy[part.df$Therapy == 'Therapied'] <- 'TX'
part.df$Therapy[part.df$Therapy == 'Relapse'] <- 'TX'
part.df$Therapy[!(part.df$Therapy %in% c('TX', 'NTX', 'Unknown'))] <- 'Unknown'
total.df$Therapy[total.df$Therapy == 'Relapse'] <- 'TX'

library(dplyr)
total_counts <- total.df %>% 
  group_by(Oncotree.Code, Therapy) %>%
  summarise(Count = n())
part_counts <- part.df %>% 
  group_by(Hugo_Symbol,Oncotree.Code, Therapy) %>%
  summarise(Count = n())

merged_counts <- inner_join(part_counts, total_counts, by = c("Oncotree.Code", "Therapy"), suffix = c("_part", "_total"))

merged_counts <- merged_counts %>%
  mutate(percentage = Count_part / Count_total*100)

library(ggplot2)
library(ggtext)

ggplot(merged_counts, aes(x = Hugo_Symbol, y = percentage, fill = interaction(Oncotree.Code, Therapy))) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = Count_part), position = position_dodge(width = 1), vjust = -0.5) +
  labs(title = "Percentage of escape genes to all samples",
       x = "Hugo Symbol",
       y = "Percentage") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values=c('gold','green','yellow','skyblue','khaki','orange','darkblue','cyan'))
ggtitle("Percentage of Count_part to Count_total")

ggplot(merged_counts, aes(x = Hugo_Symbol, y = percentage, fill = Therapy)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = Count_part), position = position_dodge(width = 0.9), vjust = -0.5, hjust = 0.5, size = 3) +
  labs(title = "Percentage of Count_part to Count_total",
       x = "Hugo Symbol",
       y = "percentage") +
  facet_wrap(.~Oncotree.Code, scales = "free_y") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5)) 
ggtitle("Percentage of Count_part to Count_total")

ggplot(total_counts, aes(x = Therapy, y = Count, fill = Oncotree.Code)) +
  geom_bar(position='fill', stat = "identity") +
  theme_bw() +
  theme(legend.position = "bottom")

ggplot(part_counts, aes(x = Therapy, y = Count, fill = Hugo_Symbol)) +
  geom_bar(position='fill', stat = "identity") +
  theme_bw() +
  #scale_fill_brewer(palette=2) +
  theme(legend.position = "bottom")

part_table <- table(part.df$Hugo_Symbol, part.df$Therapy)
result1 <- chisq.test(part_table)
print(result1)

total_table <- table(total.df$Oncotree.Code, total.df$Therapy)
result2 <- chisq.test(total_table)
print(result2)

luad1.df <- subset(total.df, Oncotree.Code=='LUAD')
luad2.df <- subset(part.df, Oncotree.Code=='LUAD')
luad2.df <- luad2.df[!duplicated(luad2.df$Tumor_Sample_Barcode),]
combined_df <- bind_rows(mutate(luad1.df, Source = "total"),mutate(luad2.df, Source = "part"))
luad_table <- table(subset(combined_df, Source == "total" | Source == "part")$Source, combined_df$Therapy)
luad_table["total", ] <- luad_table["total", ] - luad_table["part", ]##rewrite luad_total
result3 <- chisq.test(luad_table)
print(result3)

portal.df <- read.csv('d:/master project/portal samples data/bronchus and lung/outcome.csv')
portal.part <- nrow(portal.df)
portal.total <- length(list.files('d:/master project/portal samples data/bronchus and lung/'))-3-portal.part
portal1.table <- cbind(luad_table[, 1], TCGA = NA, luad_table[, 2:ncol(luad_table)])
portal1.table["part", "TCGA"] <- portal.part
portal1.table["total", "TCGA"] <- portal.total
colnames(portal1.table)[1] <- "NTX"
result4 <- chisq.test(portal1.table[,1:2])
print(result4)

portal2.table <- luad_table
portal2.table["part", "NTX"] <- luad_table["part", "NTX"] + portal.part
portal2.table['total','NTX'] <- luad_table['total','NTX'] + portal.total
result5 <- chisq.test(portal2.table[,1:2])
print(result5)
