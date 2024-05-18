library(ggplot2)
library(dplyr)
main_directory <- 'd:/master project/portal samples data/'
setwd(main_directory)
folder_names <- list.files(full.names = FALSE)
folder_names <- folder_names[!grepl('gz|10.9|outcome|antigen', folder_names, ignore.case = TRUE) & !file.exists(file.path(main_directory, "outcome"))]


merged_data <- data.frame()

counts_list <- list()

for (folder_name in folder_names) {
  folder_path <- file.path(main_directory, folder_name)
  outcome_file <- file.path(folder_path, 'outcome.csv')
  outcome <- read.csv(outcome_file)
  unique_outcome <- distinct(outcome) 
  
  num <- length(list.files(folder_path, full.names = TRUE)) - 2
  counts <- table(unique_outcome$Hugo_Symbol) / num * 100
  
  data_to_plot <- data.frame(
    Hugo_Symbol = names(counts),
    Percentage = as.numeric(counts),
    Source = folder_name
  )
  
  data_to_plot$Count <- num * data_to_plot$Percentage / 100
  
  counts_data <- data_to_plot %>%
    group_by(Hugo_Symbol, Source) %>%
    summarise(count = sum(Count)) %>%
    ungroup()
  counts_list[[folder_name]] <- counts_data
  
  merged_data <- rbind(merged_data, data_to_plot)
}

counts_data <- do.call(rbind, counts_list)

merged_data <- merge(merged_data, counts_data, by = c("Hugo_Symbol", "Source"))

ggplot(merged_data, aes(x = Hugo_Symbol, y = Percentage, fill = Source)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = count), 
            position = position_dodge(width = 0.9), 
            vjust = -0.5, size = 2) +
  labs(x = "Hugo_Symbol", y = "Percentage") +
  scale_fill_discrete(name = "Source") + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5))

