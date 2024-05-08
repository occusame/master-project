df1 <- read.csv('d:/master project/portal samples data/brain/outcome.csv')
df2 <- read.csv('d:/master project/portal samples data/breast/outcome.csv')
df3 <- read.csv('d:/master project/portal samples data/bronchus and lung/outcome.csv')
df4 <- read.csv('d:/master project/portal samples data/colon/outcome.csv')
df5 <- read.csv('d:/master project/portal samples data/hematopoietic and reticuloendothelial system/outcome.csv')
df6<- read.csv('d:/master project/portal samples data/kidney/outcome.csv')
df7 <- read.csv('d:/master project/portal samples data/ovary/outcome.csv')
df8 <- read.csv('d:/master project/portal samples data/prostate gland/outcome.csv')
library(dplyr)
library(tidyr)
merged_df <- bind_rows(df1 %>% mutate(Source = "brain"),
                       df2 %>% mutate(Source = "breast"),
                       df3 %>% mutate(Source = "bronchus\n&lung"),
                       df4 %>% mutate(Source = "colon"),
                       df5 %>% mutate(Source = "HRS"),
                       df6 %>% mutate(Source = "kidney"),
                       df7 %>% mutate(Source = "ovary"),
                       df8 %>% mutate(Source = "prostate\ngland"))
gene_counts <- merged_df %>%
  count(Source, Hugo_Symbol)
wide_gene_counts <- gene_counts %>%
  spread(key = Source, value = n, fill = 0)

print(wide_gene_counts)

library("gplots")
dt <- as.table(as.matrix(wide_gene_counts))
rownames(dt) <- paste(rownames(dt), dt[, 1], sep="")
dt <- dt[, -1, drop = FALSE]
rownames(dt) <- wide_gene_counts$Hugo_Symbol
balloonplot(t(dt), main ="mutation escape genes", xlab ="", ylab="",
            label = FALSE, show.margins = FALSE)
wide2 = wide_gene_counts
row.names(wide2) = wide2$Hugo_Symbol
wide2 = wide2[,-1]
library("graphics")
mosaicplot(wide2, shade = TRUE, las=1,
           main = "mutation escape genes")

library("vcd")
dt_numeric <- as.matrix(as.numeric(dt))
perform_chi_squared_test <- function(data) {
  test_result <- chisq.test(data)
  print(test_result)
}
perform_chi_squared_test(wide2)

