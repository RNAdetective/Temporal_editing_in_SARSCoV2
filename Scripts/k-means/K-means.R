## This script calculates overall editing levels per samples and performs K- means clustering

## This script generates figures 4A, 4B, supplementary figure 2A, 2B, table 1, supplementary (figure) tables 1, 2 and supplementary files 8, 

# libraries
library(dplyr)
library(tidyr)
library(data.table)
library(gridExtra)
library(patchwork)

setwd("H:\\Manuscript_1\\Overall_ADAR_editing_levels\\")
directoryPath <- "H:\\Manuscript_1\\Overall_ADAR_editing_levels\\"

## Overall ADAR editing levels

#Read Redi file
redi_data <- fread("Repo/TABLE1_hg19.txt", header = TRUE, fill = TRUE, stringsAsFactors = FALSE)

# Converting position column in redi_data to character (to match with coordinate in sample_file)
redi_data$Position <- as.character(redi_data$Position)

# Remove "chr" from the value in Region column of redi_data (to match with chromosome in sample_file)
redi_data$Region <- gsub("chr", "", redi_data$Region)
head(redi_data)

# Filter for reads >=100
apply_filter <- function(df) {
  df[df$stackdepth >= 100, ]
}

# calculate values for each input CSV file
calculate_values <- function(input_csv_file) {
  
  print(paste0("Processing Input CSV file : ", input_csv_file))
  sample_file <- fread(input_csv_file)
  
  #converting coordinate column in sample_file to character
  sample_file$coordinate <- as.character(sample_file$coordinate)
  print("Filtering Sample File with redi file")
  
  # merge Position and coordinate, Region and chromosome columns
  filtered_sample_file <- sample_file[.(redi_data$Region, redi_data$Position), on = c("chromosome", "coordinate"), nomatch = 0]
  print("Filtered")
  
  # filter function of >=100
  filtered_sample_file <- apply_filter(filtered_sample_file)
  
  # Convert columns G, C, and stackdepth columsn to numerical values
  filtered_sample_file <- filtered_sample_file %>%
    mutate(across(c(G, C, stackdepth), as.numeric))
  
  # This adds a new column editinglevel to filtered_sample_file 
  filtered_sample_file <- filtered_sample_file %>%
    mutate(
      editinglevel = ifelse(reference == "A", G / stackdepth, C / stackdepth)
    )
  
  filteredFileName = paste0(directoryPath,"Result\\",basename(input_csv_file),"_filtered.csv")
  write.csv(filtered_sample_file, filteredFileName, row.names = FALSE)
  
  print("Editing level added")
  
  # Calculate the sum of A, T ,
  A_count <- sum(filtered_sample_file[reference == "A",G])
  T_count <- sum(filtered_sample_file[reference == "T", C])
  coverage_sum <- sum(filtered_sample_file$stackdepth, na.rm = TRUE)
  overall_editing_level <- (A_count + T_count) / coverage_sum
  
  # Create the new data frame
  new_df <- data.frame(
    sample = basename(input_csv_file),  # Use file name as the sample name
    As = A_count,
    Ts = T_count,
    coverage = coverage_sum,
    overall_editing_level = overall_editing_level
  )
  
  print("returning new df")
  
  return(new_df)
  
}

csv_files <- list.files(paste0(directoryPath, "SourceFiles\\"), pattern = "\\.csv$", full.names = TRUE)

# Loop through each CSV file and calculate values
list_of_new_dfs <- lapply(csv_files, calculate_values)

# Concatenate all new data frames into a single data frame
new_df <- do.call(rbind, list_of_new_dfs)
View(new_df)

# Save the new_df data frame to the output file
finalResult <- file.path(directoryPath, "Result", "overall_editing.csv")
write.csv(new_df, finalResult, row.names = FALSE)

#========== ploting ============================================================

# read the final result from above that has patient_id added manually
setwd("F:\\Manuscript_1\\Overall_ADAR_editing\\Result")

Overall_ADAR <- read.csv("overall_editing.csv")
p1 <- ggplot(Overall_ADAR, aes(x = infection_status, y = overall_editing_level, color = infection_status)) +
  geom_boxplot(fill = "white", alpha = 0.5) +  
  geom_jitter(width = 0.1, alpha = 0.9) +
  labs(x = "Stages of SARS-CoV-2 infection", y = "Overall ADAR editing levels") +
  scale_color_manual(values = c("darkmagenta", "red", "blue")) +  
  theme_classic() +
  theme(axis.text.x = element_text(size = 14),  
        axis.text.y = element_text(size = 14),  
        axis.title.x = element_text(size = 10),  
        axis.title.y = element_text(size = 16),  
        legend.text = element_text(size = 14)) +
  scale_x_discrete(labels = c("Controls" = "Pre-infection", "Mid" = "Mid-infection", "Post" = "Post-infection"))

# adding stats (paired t-test)
Overall_ADAR_editing <- p1 + theme(legend.position = "none")
Overall_ADAR_editing <- Overall_ADAR_editing + stat_compare_means(comparisons = list(c("Controls", "Post"), c("Controls", "Mid")), method = "t.test", paired = TRUE)
ggsave("Overall_ADAR_editing.png", width = 8, height = 5)

#======= patientwise graph ====================================================

setwd("F:\\Manuscript_1\\Overall_ADAR_editing\\Result\\")

# result data
Overall_ADAR_patientwise <- read.csv("overall_editing _patientwise.csv", row.names = 1)
head(Overall_ADAR_patientwise)

# patient_id was added manually
patients <- rownames(Overall_ADAR_patientwise)

# calculating relative chnages  mid and post infection to pre-infection
mid_values <- Overall_ADAR_patientwise$Mid.infection - Overall_ADAR_patientwise$Pre.infection
post_values <- Overall_ADAR_patientwise$Post.infection - Overall_ADAR_patientwise$Pre.infection


# Dataframe
plot_data <- data.frame(
  Patient_ID = rep(patients, 2),
  ADAR_Editing = c(mid_values, post_values),
  Stage = rep(c("Mid", "Post"), each = length(patients))
)

# Plot
patientwisw_overall_editing <- ggplot(plot_data, aes(x = as.factor(Patient_ID), y = ADAR_Editing, fill = Stage)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  theme_minimal() +
  labs(x = "Patient ID", y = "Overall ADAR Editing Levels Mid- and Post-infection (relative to Pre-infection)",
       fill = "Infection Stage") +
  scale_fill_manual(values = c("red", "blue")) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 14),  
        axis.text.y = element_text(size = 14),  
        axis.title.x = element_text(size = 10),  
        axis.title.y = element_text(size = 16),  
        legend.text = element_text(size = 14)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black")

ggsave("patientwisw_overall_editing.png", width = 12, height = 10)


##########

library(ggplot2)
library(reshape2)

# Load data
setwd("F:\\Manuscript_1\\Overall_ADAR_editing\\Result\\")
Overall_ADAR_patientwise <- read.csv("overall_editing _patientwise.csv", row.names = 1)

# Calculating changes relative to pre-infection/baseline
baseline <- Overall_ADAR_patientwise[, 1]  # Pre-infection values
changes <- Overall_ADAR_patientwise[, -1] - baseline

# Convert changes matrix to a data frame
changes_df <- as.data.frame(t(changes))
changes_df$Patient_ID <- rownames(changes_df)
changes_df <- melt(changes_df, id.vars = "Patient_ID")
colnames(changes_df) <- c("Infection_stage", "Patient_ID", "Overall_ADAR_editing_levels")

# Create ggplot
Overall_ADAR <- ggplot(changes_df, aes(x = Patient_ID, y = Overall_ADAR_editing_levels, fill = Infection_stage)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Patient_ID", y = "Overall ADAR editing levels - relative to pre-infection") +
  theme_classic() +
  theme(legend.position = "bottom") +
  theme(axis.text.x = element_text(size = 12),  
        axis.text.y = element_text(size = 12),  
        axis.title.x = element_text(size = 16),  
        axis.title.y = element_text(size = 16),  
        legend.text = element_text(size = 14))  

Overall_ADAR <- ggplot(changes_df, aes(x = Patient_ID, y = Overall_ADAR_editing_levels, fill = Infection_stage)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Patient_ID", y = "Overall ADAR editing levels - relative to pre-infection") +
  theme_classic() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(size = 12,angle = 45, vjust = 1, hjust = 1),  # Adjust angle here
        axis.text.y = element_text(size = 12),  
        axis.title.x = element_text(size = 16),  
        axis.title.y = element_text(size = 16),  
        legend.text = element_text(size = 14))



# Create ggplot
Overall_ADAR <- ggplot(changes_df, aes(x = Overall_ADAR_editing_levels, y = Patient_ID, fill = Infection_stage)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Overall ADAR editing levels - relative to pre-infection", y = "Patient_ID") +
  theme_classic() +
  theme(legend.position = "bottom") +
  theme(axis.text.x = element_text(size = 12),  
        axis.text.y = element_text(size = 12),  
        axis.title.x = element_text(size = 14),  
        axis.title.y = element_text(size = 14),  
        legend.text = element_text(size = 12)) 


#========================  K-means =============================================

# library
library(ggplot2) 
library(ggfortify) #(#https://www.rdocumentation.org/packages/ggbio/versions/1.20.1/topics/autoplot)
library(dplyr)

# Pre vs Mid
setwd("F:\\Manuscript_1\\Overall_ADAR_editing\\K_means")

# read data
data <- read.csv("K_means_data_pre_post.csv", header = TRUE, sep = ",")
head(data)

# data labels
data_labels_patientID <- data$patient_id

# k means data
data_K_means <- select(data, c(2,3))
head(data_K_means)

# plotting to get overview of data 
plot(data_K_means, pch=1, cex = 2) 

# WSS - within group sum of squares
wssplot <- function(data, nc=15, seed=1234){
  wss <- (nrow(data)-1)*sum(apply(data,2,var))
  for (i in 2:nc){
    set.seed(seed)
    wss[i] <- sum(kmeans(data, centers=i)$withinss)}
  plot(1:nc, wss, type="b", xlab="Number of Clusters",
       ylab="Within groups sum of squares")
  wss
}
wssplot(data_K_means)

# k-means cluster
KM <- kmeans(data_K_means, 2, nstart = 1000)
KM

# Assigned clusters 
Km.clusters <- KM$cluster
Km.clusters

# Adding patient Id to the data
row.names(data_K_means) <- paste(data_labels_patientID, 1:dim(data_K_means)[1])
head(data_K_means)

# Plotting clustering result
plot(data_K_means, col = KM$cluster, pch = 20, cex = 3)     

# Find centers
KM$centers


## extract cluster data
# cluster1
C1 <- which(KM$cluster == 1)
Cluster.1 <- data[C1,]
#dim(Cluster.1)
write.csv(Cluster.1, "post_cluster1.csv")

#cluster2
C2 <- which(KM$cluster == 2)
Cluster.2 <- data[C2,]
#dim(Cluster.2)
write.csv(Cluster.2, "post_cluster2.csv")

