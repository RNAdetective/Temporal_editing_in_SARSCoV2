## This script compares expression of ADARs and ADAR isoforms across stage of SARS-CoV-2 infection. 

## This generates figures 1A, 1B, 1C, 2A, and 2B, and supplementary tables - 3A, 3B, 3C, 4A, and 4B

library(dplyr)
library(stringr)
library(ggpubr)
library(ggplot2)
library(gridExtra)


# Expression of ADAR enzymes
setwd ("F:\\Manuscript_1\\Counts\\Counts")

#================== Pre-infection ===========================

# Import pre-infection files
file_list <- list.files(pattern = "_c.tab$")
data_list <- list()
for (file in file_list) {
  data <- read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    data_list[[file]] <- data
}

# Combine all the above pre-infection count files to a master dataframe
Pre_infection_df <- do.call(rbind, data_list)
head(Pre_infection_df)

# Add column for infection stage 
Pre_infection_df <- cbind(Pre_infection_df, "Infection_stage" = "Pre-infection")

# Filter ADAR genes 
Pre_ADAR1 <- subset(Pre_infection_df, Gene.Name == 'ADAR')
Pre_ADAR2 <- subset(Pre_infection_df, Gene.Name == 'ADARB1')
Pre_ADAR3 <- subset(Pre_infection_df, Gene.Name == 'ADARB2')

#================== Mid-infection =========================

# Import Mid-infection files
file_list <- list.files(pattern = "_m.tab$")
data_list <- list()
for (file in file_list) {
  data <- read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  data_list[[file]] <- data
}

# combine all the above Mid-infection count files to a master dataframe
Mid_infection_df <- do.call(rbind, data_list)
head(Mid_infection_df)

# Add column for infection stage 
Mid_infection_df <- cbind(Mid_infection_df, "Infection_stage" = "Mid-infection")

# filter ADAR genes 
Mid_ADAR1 <- subset(Mid_infection_df, Gene.Name == 'ADAR')
Mid_ADAR2 <- subset(Mid_infection_df, Gene.Name == 'ADARB1')
Mid_ADAR3 <- subset(Mid_infection_df, Gene.Name == 'ADARB2')

#================ Post-infection=================

# Import Post-infection files
file_list <- list.files(pattern = "_p.tab$")
data_list <- list()
for (file in file_list) {
  data <- read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  data_list[[file]] <- data
}

# combine all the above Post-infection count files to a master dataframe
Post_infection_df <- do.call(rbind, data_list)
head(Post_infection_df)

# Add column for infection stage 
Post_infection_df <- cbind(Post_infection_df, "Infection_stage" = "Post-infection")

# filter ADAR genes 
Post_ADAR1 <- subset(Post_infection_df, Gene.Name == 'ADAR')
Post_ADAR2 <- subset(Post_infection_df, Gene.Name == 'ADARB1')
Post_ADAR3 <- subset(Post_infection_df, Gene.Name == 'ADARB2')

# combine pre,mid, and post ADAR1
ADAR1_df <- rbind(Pre_ADAR1, Mid_ADAR1, Post_ADAR1)
head(ADAR1_df)
# combine pre,mid, and post ADAR2
ADAR2_df <- rbind(Pre_ADAR2, Mid_ADAR2, Post_ADAR2)

# combine pre,mid, and post ADAR3
ADAR3_df <- rbind(Pre_ADAR3, Mid_ADAR3, Post_ADAR3)


## ADAR1 plot

# Reorder infection stage 
ADAR1_df$Infection_stage <- factor(ADAR1_df$Infection_stage, levels = c("Pre-infection", "Mid-infection", "Post-infection"))
p1 <- ggplot(ADAR1_df, aes(x = Infection_stage, y = TPM, color = Infection_stage)) +
  geom_boxplot(fill = "white", alpha = 0.5) +  
  geom_jitter(width = 0.1, alpha = 0.9) +
  labs(x = "Stages of SARS-CoV-2 infection", y = "ADAR1 expression in TPM") +
  scale_color_manual(values = c("darkmagenta", "red", "blue")) +  
  theme_classic() +
  theme(axis.text.x = element_text(size = 14),  
        axis.text.y = element_text(size = 14),  
        axis.title.x = element_text(size = 10),  
        axis.title.y = element_text(size = 16),  
        legend.text = element_text(size = 8))

# add stats (paired t-test)
ADAR1 <- p1 + theme(legend.position = "none")
ADAR1 <- ADAR1 + stat_compare_means(comparisons = list(c("Pre-infection", "Post-infection"), c("Pre-infection", "Mid-infection")), method = "t.test", paired = TRUE)


## ADAR2 plot

# Reorder infection stage 
ADAR2_df$Infection_stage <- factor(ADAR2_df$Infection_stage, levels = c("Pre-infection", "Mid-infection", "Post-infection"))

p2 <- ggplot(ADAR2_df, aes(x = Infection_stage, y = TPM, color = Infection_stage)) +
  geom_boxplot(fill = "white", alpha = 0.5) +  
  geom_jitter(width = 0.1, alpha = 0.9) +
  labs(x = "Stages of SARS-CoV-2 infection", y = "ADAR2 expression in TPM") +
  scale_color_manual(values = c("darkmagenta", "red", "blue")) +  
  theme_classic() +
  theme(axis.text.x = element_text(size = 14),  
        axis.text.y = element_text(size = 14),  
        axis.title.x = element_text(size = 10),  
        axis.title.y = element_text(size = 16),  
        legend.text = element_text(size = 14))

# add stats (paired t-test)
ADAR2 <- p2 + theme(legend.position = "none")
ADAR2 <- ADAR2 + stat_compare_means(comparisons = list(c("Pre-infection", "Post-infection"), c("Pre-infection", "Mid-infection")), method = "t.test", paired = TRUE)



## ADAR3 plot

# Reorder infection stage 
ADAR3_df$Infection_stage <- factor(ADAR3_df$Infection_stage, levels = c("Pre-infection", "Mid-infection", "Post-infection"))

p3 <- ggplot(ADAR3_df, aes(x = Infection_stage, y = TPM, color = Infection_stage)) +
  geom_boxplot(fill = "white", alpha = 0.5) +  
  geom_jitter(width = 0.1, alpha = 0.9) +
  labs(x = "Stages of SARS-CoV-2 infection", y = "ADAR3 expression in TPM") +
  scale_color_manual(values = c("darkmagenta", "red", "blue")) +  
  theme_classic() +
  theme(axis.text.x = element_text(size = 14),  
        axis.text.y = element_text(size = 14),  
        axis.title.x = element_text(size = 10),  
        axis.title.y = element_text(size = 16),  
        legend.text = element_text(size = 14))
  


# add stats (paired t-test)
ADAR3 <- p3 + theme(legend.position = "none")
ADAR3 <- ADAR3 + stat_compare_means(comparisons = list(c("Pre-infection", "Post-infection"), c("Pre-infection", "Mid-infection")), method = "t.test", paired = TRUE)

# Arrange and label in grid
ggarrange(ADAR1, ADAR2, ADAR3, ncol = 3, labels= c("A", "B", "C"))


