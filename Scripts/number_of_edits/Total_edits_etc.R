## This script calculates the total number of ADAR edits per samples, unique edits per stage of infection, regionwise distribution of unique edits, and reactome pathways (results section alterations in number and genomic distribution of ADAR editing sites)

## This generates figures 3A, 3B, 3C, and 3D, supplementary figures 1A, 1B, and supplementary tables 5, 6A,6B, 6C, 6D, 6E, and 7

##libraries
library(dplyr)
library(stringr)
library(ggplot2)
library(ggpubr) 
library(reshape2)
library(cowplot)
library(gridExtra)
#install.packages("UpSetR")
library(UpSetR)
#BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
library(patchwork)
library(RColorBrewer)
library(reshape2)


# Directory
setwd("F:\\Manuscript_1\\Putative_edits_ATGC")

#------ Total number of ADAR (A>G and T>C edits per sample-------------------------------------------------------


# This function will run on each file "*filtered_snps_finalAll.csv" and apply the filters
process_file <- function(file_path) {
  
  VCF_files <- read.csv(file_path, header = TRUE, sep = ",")
  
  # Filter1: keep rows with "." in "ID" column, this is to remove the rows with SNP ids
  filtered_data <-  VCF_files %>%
    filter(ID == ".")
  
  # Filter2: remove rows that have 2 values in column "ALT" - ie these are polymorphic sites
  filtered_data <- filtered_data %>%
    filter(!str_detect(ALT, ","))
  
  # Filter3: keep rows where either "REF = A" and "ALT = G" or "REF = T" and "ALT = C" - to keep only ADAR edited sites
  filtered_data <- filtered_data %>%
    filter((REF == "A" & ALT == "G") | (REF == "T" & ALT == "C"))
  
  # Count total number of A and Ts - ie counting number of ADAR edits
  ref_tally_A <- sum(filtered_data$REF == "A") # this will include all the A to G
  ref_tally_T <- sum(filtered_data$REF == "T") # will include all the T to c
  
  # Extract the file name without extension to identify files in the final dataframe (manually check to re-confirm)
  file_name <- tools::file_path_sans_ext(basename(file_path)) ## tools package
  
  return(list(ref_tally_A = ref_tally_A, ref_tally_T = ref_tally_T, file_name = file_name))
}

# Directory containing the files
directory_path <- "F:\\Manuscript_1\\Putative_edits_ATGC"

#empty list for the results 
results_list <- list()

# loop 
files <- list.files(pattern = "*.csv", full.names = TRUE, path = directory_path)
for (file_path in files) {
  result <- process_file(file_path)
  results_list[[length(results_list) + 1]] <- result # to ensure that the result is added to the end of the list
}

# Extracting vectors of A T and file names 
ref_tally_A_vector <- sapply(results_list, function(result) result$ref_tally_A)
ref_tally_T_vector <- sapply(results_list, function(result) result$ref_tally_T)
file_names_vector <- sapply(results_list, function(result) result$file_name)

# Join all the above vectors of number of A and t edits to a dataframe
Total_AG_TC <- data.frame(sample = file_names_vector, Total_AG = ref_tally_A_vector, Total_TC = ref_tally_T_vector, 
                          TOTAL_ADAR = ref_tally_A_vector + ref_tally_T_vector)
head(Total_AG_TC)

# write the above adtaframe
write.csv(Total_AG_TC, "Result_Total_number_edits(ATGC).csv", row.names = FALSE)


# ------------------ Box plot comparing total number of putative ADAR edits within the cohort-----------------------

Total_edits <- read.csv("Result_Total_number_edits(ATGC).csv")
# infection status was added manually

Total_edits_plot <- ggplot(Total_edits, aes(x = infection.status, y = TOTAL_ADAR, color = infection.status)) +
  geom_boxplot(fill = "white", alpha = 0.5) +  
  geom_jitter(width = 0.1, alpha = 0.9) +
  labs(x = "Stages of SARS-CoV-2 infection", y = "Total number of putative ADAR edits") +
  scale_color_manual(values = c("darkmagenta", "red", "blue")) +  
  theme_classic() +
  theme(axis.text.x = element_text(size = 14),  
        axis.text.y = element_text(size = 14),  
        axis.title.x = element_text(size = 10),  
        axis.title.y = element_text(size = 16),  
        legend.text = element_text(size = 14)) +
  scale_x_discrete(labels = c("Controls" = "Pre-infection", "Mid" = "Mid-infection", "Post" = "Post-infection"))

# adding stats (paired t-test)
Total_edits_plot <- Total_edits_plot + theme(legend.position = "none")
Total_edits_plot <- Total_edits_plot + stat_compare_means(comparisons = list(c("Controls", "Post"), c("Controls", "Mid")), method = "t.test", paired = TRUE)

ggsave("Total_edits_plot.png", width = 8, height = 5)

#=================  patientwise changes in total number of edits ======================

# make Total_edits long format
Total_edits <- Total_edits %>%
  select(Patient_Id, infection.status, TOTAL_ADAR) %>%
  spread(key = infection.status, value = TOTAL_ADAR)

write.csv(Total_edits, "F:\\Manuscript_1\\Putative_edits_ATGC\\patientwise_putativeedits.csv")

Total_edits <- read.csv("Result_Total_number_edits(ATGC).csv")
average_matrix <- read.csv ("F:\\Manuscript_1\\Alignment_matrix.csv", sep = "\t")

#add total nunmber of aligend reads to the "total edits dataframe"

# Check if the sample column values are same/order
if (identical(Total_edits$sample, average_matrix$Sample_SRR)) {
  
  Total_edits$aligned_reads <- average_matrix$Number.of.aligned.reads
  head(Total_edits)
  
  # save df with aligned reads
  write.csv(Total_edits, "Total_edits_with_aligned_reads.csv", row.names = FALSE)
  
} else {
  print("The sample IDs are diff")
}

Total_edits_with_aligned_reads <- read.csv("Total_edits_with_aligned_reads.csv")

# normalized ADAR editing per sample
Total_edits_with_aligned_reads$normalized_ADAR <- (Total_edits_with_aligned_reads$TOTAL_ADAR / Total_edits_with_aligned_reads$aligned_reads) * 1000000
head(Total_edits_with_aligned_reads)
write.csv(Total_edits_with_aligned_reads, "Normalized_Total_edits.csv", row.names = FALSE)

# Filter data as per infection stage 
baseline <- Total_edits_with_aligned_reads %>%
  filter(infection.status == "Controls") %>%
  select(Patient_Id, normalized_ADAR) %>%
  rename(Pre_ADAR = normalized_ADAR)  # Rename column for clarity
Mid <- Total_edits_with_aligned_reads %>%
  filter(infection.status == "Mid") %>%
  select(Patient_Id, normalized_ADAR) %>%
  rename(Mid_ADAR = normalized_ADAR)  # Rename column for clarity

Post <- Total_edits_with_aligned_reads %>%
  filter(infection.status == "Post") %>%
  select(Patient_Id, normalized_ADAR) %>%
  rename(Post_ADAR = normalized_ADAR)  # Rename column for clarity

# Merge Pre, Mid, and Post based on Patient_Id
merged_df <- baseline %>%
  left_join(Mid, by = "Patient_Id") %>%
  left_join(Post, by = "Patient_Id")

# Reshape the df
plot_data <- merged_df %>%
  select(Patient_Id, Pre_ADAR, Mid_ADAR, Post_ADAR) %>%
  melt(id.vars = "Patient_Id", variable.name = "Stage", value.name = "ADAR")

# Convert Stage to factor to ensure correct order
plot_data$Stage <- factor(plot_data$Stage, levels = c("Pre_ADAR", "Mid_ADAR", "Post_ADAR"))

# This will calculate changes releatuve to pre-infection
relative_changes <- merged_df %>%
  mutate(
    Mid_Change = Mid_ADAR - Pre_ADAR,
    Post_Change = Post_ADAR - Pre_ADAR
  ) %>%
  select(Patient_Id, Mid_Change, Post_Change)

# Reshape the data for plot
plot_data <- melt(relative_changes, id.vars = "Patient_Id", 
                  variable.name = "Stage", value.name = "Absolute_Change")


#plot- patientwise
Total_edits_patientwise_normalized <- ggplot(plot_data, aes(x = as.factor(Patient_Id), y = Absolute_Change, fill = Stage)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  theme_minimal() +
  labs(x = "Patient ID",
       y = "Total number of putative ADAR edits - relative to pre-infection",
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


ggsave("Total_edits_patientwise_normalized.png", width = 13, height = 9)

#========================  sites consistently edited across stages of infection =========

#working directory
directory <- "F:\\Manuscript_1\\Frequency_Counts\\Frequency_Counts"
setwd(directory)

# list if all frequency counts files
csv_files <- dir(pattern = "\\.csv$")

# Loop through each CSV file
for (file in csv_files) {
  data <- read.csv(file, stringsAsFactors = FALSE)
  
  filtered_data <- data %>%
    filter(reference == "A" | reference == "T")
  
  filtered_data <- filtered_data %>%
    mutate(editing_level = ifelse(reference == "A", G / stackdepth, C / stackdepth))
  
  filtered_data <- filtered_data %>%
    filter(stackdepth >= 100)
  
  write.csv(filtered_data, file.path("filtered_files", paste0("filtered_", file)), row.names = FALSE)
}

################################################################################
################################################################################

# Filter "consistently edited sites" 
directory <- "F:\\Manuscript_1\\Frequency_Counts\\Frequency_Counts\\filtered_files"
setwd(directory)

## sites shared across 20% of pre infection samples, stack depth >=100
file_names <- c("filtered_SRR18302396All.csv",
                "filtered_SRR18302665All.csv",
                "filtered_SRR18302708All.csv",
                "filtered_SRR18302720All.csv",
                "filtered_SRR18302738All.csv",
                "filtered_SRR18302745All.csv",
                "filtered_SRR18302752All.csv",
                "filtered_SRR18302764All.csv",
                "filtered_SRR18302780All.csv",
                "filtered_SRR18302865All.csv",
                "filtered_SRR18302873All.csv",
                "filtered_SRR18303008All.csv",
                "filtered_SRR18303080All.csv",
                "filtered_SRR18303098All.csv",
                "filtered_SRR18303106All.csv",
                "filtered_SRR18303112All.csv",
                "filtered_SRR18303118All.csv",
                "filtered_SRR18303180All.csv",
                "filtered_SRR18303189All.csv",
                "filtered_SRR18303195All.csv",
                "filtered_SRR18303427All.csv",
                "filtered_SRR18303450All.csv",
                "filtered_SRR18303457All.csv",
                "filtered_SRR18303561All.csv",
                "filtered_SRR18303602All.csv",
                "filtered_SRR18303610All.csv",
                "filtered_SRR18303617All.csv",
                "filtered_SRR18303633All.csv",
                "filtered_SRR18303640All.csv",
                "filtered_SRR18303646All.csv",
                "filtered_SRR18303656All.csv",
                "filtered_SRR18303662All.csv",
                "filtered_SRR18303668All.csv",
                "filtered_SRR18303680All.csv",
                "filtered_SRR18303686All.csv",
                "filtered_SRR18303718All.csv",
                "filtered_SRR18303724All.csv",
                "filtered_SRR18303734All.csv",
                "filtered_SRR18303761All.csv",
                "filtered_SRR18303785All.csv",
                "filtered_SRR18303834All.csv",
                "filtered_SRR18303915All.csv",
                "filtered_SRR18303925All.csv",
                "filtered_SRR18303970All.csv",
                "filtered_SRR18304012All.csv")


# store the frequency files
dataframes_list <- list()

# This Loop will read and filter each file
for (file_name in file_names) {
  # Read the files
  VCF_data <- read.csv(file_name, header = TRUE, sep = ",", stringsAsFactors = FALSE)
  
  # Step 1: Concatenate 'chromosome' and 'coordinate' columns with "_"
  VCF_data$CHROM_POS <- paste(VCF_data$chromosome, VCF_data$coordinate, sep="_")
  
  # Step 2: Filter rows with 'stackdepth' greater than or equal to 100
  VCF_data <- VCF_data %>%
    filter(stackdepth >= 100)
  
  # Step 3: Filter rows with values "A" or "T" in the column "reference"
  VCF_data <- VCF_data %>%
    filter(reference == 'A' | reference == 'T')
  
  # Add the filtered data frame to the list
  dataframes_list[[file_name]] <- VCF_data
}
View(head(VCF_data))


# This will find common values present in at least 20% of input files
min_presence <- 0.2 * length(dataframes_list) # change the value here for a different percentage
common_values_all <- names(Filter(function(x) x >= min_presence, table(unlist(lapply(dataframes_list, function(df) df$CHROM_POS)))))
common_values_all # the way this is calculated is it will ensure that each "CHROM_POS" pair is present in at least 20% of samples

# Write the output 
write.csv(data.frame(merged_column = common_values_all), file = "Consistent_edits_pre_20%samples.csv", row.names = FALSE)

#============================ mid-infection ======================================================================================

## sites shared across 20% of Mid-infection samples, stack depth >=100
file_names <- c("filtered_SRR18302393All.csv",
                "filtered_SRR18302423All.csv",
                "filtered_SRR18302666All.csv",
                "filtered_SRR18302709All.csv",
                "filtered_SRR18302721All.csv",
                "filtered_SRR18302739All.csv",
                "filtered_SRR18302746All.csv",
                "filtered_SRR18302753All.csv",
                "filtered_SRR18302765All.csv",
                "filtered_SRR18302781All.csv",
                "filtered_SRR18302808All.csv",
                "filtered_SRR18302948All.csv",
                "filtered_SRR18303011All.csv",
                "filtered_SRR18303081All.csv",
                "filtered_SRR18303100All.csv",
                "filtered_SRR18303107All.csv",
                "filtered_SRR18303113All.csv",
                "filtered_SRR18303119All.csv",
                "filtered_SRR18303181All.csv",
                "filtered_SRR18303190All.csv",
                "filtered_SRR18303196All.csv",
                "filtered_SRR18303428All.csv",
                "filtered_SRR18303458All.csv",
                "filtered_SRR18303562All.csv",
                "filtered_SRR18303603All.csv",
                "filtered_SRR18303611All.csv",
                "filtered_SRR18303618All.csv",
                "filtered_SRR18303634All.csv",
                "filtered_SRR18303641All.csv",
                "filtered_SRR18303647All.csv",
                "filtered_SRR18303657All.csv",
                "filtered_SRR18303663All.csv",
                "filtered_SRR18303669All.csv",
                "filtered_SRR18303681All.csv",
                "filtered_SRR18303687All.csv",
                "filtered_SRR18303719All.csv",
                "filtered_SRR18303727All.csv",
                "filtered_SRR18303735All.csv",
                "filtered_SRR18303801All.csv",
                "filtered_SRR18303813All.csv",
                "filtered_SRR18303819All.csv",
                "filtered_SRR18303866All.csv",
                "filtered_SRR18303871All.csv",
                "filtered_SRR18303926All.csv",
                "filtered_SRR18304013All.csv")

#An empty list is opened to store dataframes
dataframes_list <- list()

## This Loop will read and filter each file
for (file_name in file_names) {
  # Read the files
  VCF_data <- read.csv(file_name, header = TRUE, sep = ",", stringsAsFactors = FALSE)
  
  # Step 1: Concatenate 'chromosome' and 'coordinate' columns with "_"
  VCF_data$CHROM_POS <- paste(VCF_data$chromosome, VCF_data$coordinate, sep="_")
  
  # Step 2: Filter rows with 'stackdepth' greater than or equal to 100
  VCF_data <- VCF_data %>%
    filter(stackdepth >= 100)
  
  # Step 3: Filter rows with values "A" or "T" in the column "reference"
  VCF_data <- VCF_data %>%
    filter(reference == 'A' | reference == 'T')
  
  # Add the filtered data frame to the list
  dataframes_list[[file_name]] <- VCF_data
}
View(head(VCF_data))

# This will find common values present in at least 20% of input files
min_presence <- 0.2 * length(dataframes_list) # change the value here for a different percentage
common_values_all <- names(Filter(function(x) x >= min_presence, table(unlist(lapply(dataframes_list, function(df) df$CHROM_POS)))))
common_values_all # the way this is calculated is it will ensure that each "CHROM_POS" pair is present in at least 80% of samples

# Write 
write.csv(data.frame(merged_column = common_values_all), file = "Consistent_edits_mid_20%samples.csv", row.names = FALSE)

#======================= post infection ================================================================================================================

# sites shared across 20% of Post-infection samples, stack depth >=100
file_names <- c("filtered_SRR18302673All.csv",
                "filtered_SRR18302711All.csv",
                "filtered_SRR18302725All.csv",
                "filtered_SRR18302743All.csv",
                "filtered_SRR18302750All.csv",
                "filtered_SRR18302757All.csv",
                "filtered_SRR18302767All.csv",
                "filtered_SRR18302783All.csv",
                "filtered_SRR18302934All.csv",
                "filtered_SRR18302941All.csv",
                "filtered_SRR18302971All.csv",
                "filtered_SRR18303015All.csv",
                "filtered_SRR18303086All.csv",
                "filtered_SRR18303105All.csv",
                "filtered_SRR18303111All.csv",
                "filtered_SRR18303117All.csv",
                "filtered_SRR18303124All.csv",
                "filtered_SRR18303184All.csv",
                "filtered_SRR18303194All.csv",
                "filtered_SRR18303198All.csv",
                "filtered_SRR18303223All.csv",
                "filtered_SRR18303266All.csv",
                "filtered_SRR18303280All.csv",
                "filtered_SRR18303431All.csv",
                "filtered_SRR18303453All.csv",
                "filtered_SRR18303463All.csv",
                "filtered_SRR18303567All.csv",
                "filtered_SRR18303608All.csv",
                "filtered_SRR18303615All.csv",
                "filtered_SRR18303622All.csv",
                "filtered_SRR18303639All.csv",
                "filtered_SRR18303645All.csv",
                "filtered_SRR18303652All.csv",
                "filtered_SRR18303661All.csv",
                "filtered_SRR18303667All.csv",
                "filtered_SRR18303674All.csv",
                "filtered_SRR18303685All.csv",
                "filtered_SRR18303691All.csv",
                "filtered_SRR18303721All.csv",
                "filtered_SRR18303729All.csv",
                "filtered_SRR18303736All.csv",
                "filtered_SRR18303917All.csv",
                "filtered_SRR18303930All.csv",
                "filtered_SRR18303972All.csv",
                "filtered_SRR18304019All.csv")


#An empty list is opened to store dataframes
dataframes_list <- list()

## This Loop will read and filter each file
for (file_name in file_names) {
  # Read the files
  VCF_data <- read.csv(file_name, header = TRUE, sep = ",", stringsAsFactors = FALSE)
  
  # Step 1: Concatenate 'chromosome' and 'coordinate' columns with "_"
  VCF_data$CHROM_POS <- paste(VCF_data$chromosome, VCF_data$coordinate, sep="_")
  
  # Step 2: Filter rows with 'stackdepth' greater than or equal to 100
  VCF_data <- VCF_data %>%
    filter(stackdepth >= 100)
  
  # Step 3: Filter rows with values "A" or "T" in the column "reference"
  VCF_data <- VCF_data %>%
    filter(reference == 'A' | reference == 'T')
  
  # Add the filtered data frame to the list
  dataframes_list[[file_name]] <- VCF_data
}


# Find common values present in at least 20% of input files
min_presence <- 0.2 * length(dataframes_list) # change the value here for a different percentage
common_values_all <- names(Filter(function(x) x >= min_presence, table(unlist(lapply(dataframes_list, function(df) df$CHROM_POS)))))
common_values_all # the way this is calculated is it will ensure that each "CHROM_POS" pair is present in at least 80% of samples

# Write 
write.csv(data.frame(merged_column = common_values_all), file = "Consistent_edits_post_20%samples.csv", row.names = FALSE)

#======================================== unique an shared sites =====================================================================================

setwd("F:\\Manuscript_1\\Frequency_Counts\\Frequency_Counts\\filtered_files - Copy")
# Read the data
pre_infection <- read.csv("Consistent_edits_pre_20%samples.csv")
mid_infection <- read.csv("Consistent_edits_mid_20%samples.csv")
post_infection <- read.csv("Consistent_edits_post_20%samples.csv")

# Create a list containing the editing sites from each stage
#editing_sites_list <- list(
#pre = pre_infection$merged_column,
# mid = mid_infection$merged_column,
# post = post_infection$merged_column
#)

# Create the upset plot
#upset(fromList(editing_sites_list), 
# nintersects = 40,
# nsets = 6, 
#  order.by = "freq", 
# main.bar.color = "skyblue", sets.bar.color = "skyblue", 
# text.scale = c(1.5, 1.2),
# point.size = 3.5, 
# line.size = 1)


# ============= unique and shared consistently edited sites =========================================================================
# code adopted from https://github.com/cxli233/customized_upset_plots
my_list <- list(
  pre_infection = pre_infection$merged_column, 
  Mid_infection = mid_infection$merged_column, 
  Post_infection = post_infection$merged_column)

comb_mat <- make_comb_mat(my_list)
my_names <- set_name(comb_mat)

my_set_sizes <- set_size(comb_mat) %>% 
  as.data.frame() %>% 
  rename(sizes = ".") %>% 
  mutate(Set = row.names(.)) 

p1 <- my_set_sizes %>% 
  mutate(Set = reorder(Set, sizes)) %>% 
  ggplot(aes(x = Set, y= sizes)) +
  geom_bar(stat = "identity", aes(fill = Set), alpha = 0.8, width = 0.7) +
  geom_text(aes(label = sizes), 
            size = 5, angle = 90, hjust = 0, y = 1) +
  #scale_fill_manual(values = brewer.pal(4, "Set2"),  # feel free to use some other colors  
  scale_fill_manual(values = c("darkmagenta", "red", "blue"),                 
                    limits = my_names) + 
  labs(x = NULL,
       y = "Set size",
       fill = NULL) +
  theme_classic() +
  theme(legend.position = "right",
        text = element_text(size= 14),
        axis.ticks.y = element_blank(),
        axis.text = element_blank()
  )


get_legend <- function(p) {
  tmp <- ggplot_gtable(ggplot_build(p))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}

p2 <- get_legend(p1)


my_overlap_sizes <- comb_size(comb_mat) %>% 
  as.data.frame() %>% 
  rename(overlap_sizes = ".") %>% 
  mutate(category = row.names(.))

p3 <- my_overlap_sizes %>% 
  mutate(category = reorder(category, -overlap_sizes)) %>% 
  ggplot(aes(x = category, y = overlap_sizes)) +
  geom_bar(stat = "identity", fill = "darkturquoise", color = NA, alpha = 0.8, width = 0.7) +
  geom_text(aes(label = overlap_sizes, y = 0), 
            size = 5, hjust = 0, vjust = 0.5) +
  labs(y = "Intersect sizes",
       x = NULL) +
  theme_classic() +
  theme(text = element_text(size= 14, color = "black"),
        axis.text =element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_text(hjust = 0),
  ) +
  coord_flip()


my_overlap_matrix <- str_split(string = my_overlap_sizes$category, pattern = "", simplify = T) %>% 
  as.data.frame() 

colnames(my_overlap_matrix) <- my_names

my_overlap_matrix_tidy <- my_overlap_matrix %>% 
  cbind(category = my_overlap_sizes$category) %>% 
  pivot_longer(cols = !category, names_to = "Set", values_to = "value") %>% 
  full_join(my_overlap_sizes, by = "category") %>% 
  full_join(my_set_sizes, by = "Set")

p4 <- my_overlap_matrix_tidy %>% 
  mutate(category = reorder(category, -overlap_sizes)) %>%  
  mutate(Set = reorder(Set, sizes)) %>%  
  ggplot(aes(x = Set, y = category))+
  geom_tile(aes(fill = Set, alpha = value), color = "grey30", size = 1) +
  scale_fill_manual(values = c("darkmagenta", "red", "blue"), 
                    limits = my_names) +
  scale_alpha_manual(values = c(0.8, 0),  # color the grid for 1, don't color for 0. 
                     limits = c("1", "0")) +
  labs(x = "Sets",  
       y = "Overlap") +
  theme_minimal() +
  theme(legend.position = "none",
        text = element_text(color = "black", size= 14),
        panel.grid = element_blank(),
        axis.text = element_blank()
  )

unique_shared <- wrap_plots(p1, p2, p4, p3, 
                            nrow = 2, 
                            ncol = 2,
                            heights = c(1, 2), # the more rows in the lower part, the longer it should be
                            widths = c(1, 0.8),
                            guides = "collect") &
  theme(legend.position = "none") 

unique_shared <- unique_shared + plot_annotation(title = "B")
ggsave("unique_shared.png", width = 8, height = 5)


#========================================================================================

# sites unique to the three stages of infection were filterd using an online venn diagram maker 


#==================== annotating with rediportal =========================================

redi <- read.delim("TABLE1_hg19.txt", header = TRUE) # REDI portal ADAR editing sites 
head(redi)

# remove 'chr' an combine position and region in redi
redi$Region <- sub("chr", "", redi$Region)

redi$Combined_Coordinate <- paste(redi$Region, redi$Position, sep = "_")

Pre_unique <- read.csv("Pre_unique.csv", header = TRUE)
Mid_unique <- read.csv("Mid_unique.csv", header = TRUE)
Post_unique <- read.csv("Post_unique.csv", header = TRUE)

# ==pre
Pre_unique <- redi[redi$Combined_Coordinate %in% Pre_unique$Pre_unique, ]
View(Pre_unique)
write.csv(Pre_unique, file = "Pre_unique_annotated.csv", col.names = TRUE)

# == mid
Mid_unique <- redi[redi$Combined_Coordinate %in% Mid_unique$Mid_unique, ]
View(Mid_unique)
write.csv(Mid_unique, file = "Mid_unique_annotated.csv", col.names = TRUE)

# ==post

Post_unique <- redi[redi$Combined_Coordinate %in% Post_unique$Post_unique, ]
View(Post_unique)
write.csv(Post_unique, file = "Post_unique_annotated.csv", col.names = TRUE)



#=========== genic region====================================================================================================

Pre <- read.csv("Pre_unique_annotated.csv", header = TRUE)
Pre$infection_stage <- 'Pre-infection' # adding column
dim(Pre)
colnames(Pre)
Mid <- read.csv("Mid_unique_annotated.csv", header = TRUE)
dim(Mid)
Mid$infection_stage <- 'Mid-infection' # adding column
colnames(Mid)

Post <- read.csv("Post_unique_annotated.csv", header = TRUE)
Post$infection_stage <- 'Post-infection' # adding column
dim(Post)


combined_df <- rbind(Pre, Mid, Post)
View(combined_df)


genic_counts <- combined_df %>%
  filter(!is.na(AAChange.wgEncodeGencodeBasicV34lift37)) %>%
  group_by(infection_stage, AAChange.wgEncodeGencodeBasicV34lift37) %>%
  summarise(counts = n()) %>%
  ungroup()

# Reorder the levels 
genic_counts$infection_stage <- factor(genic_counts$infection_stage, levels = c("Pre-infection", "Mid-infection", "Post-infection"))

# Plotting 
genic_region <- ggplot(genic_counts, aes(x = counts, y = AAChange.wgEncodeGencodeBasicV34lift37, fill = infection_stage)) +
  geom_bar(stat = "identity", position = "dodge") +  # Bars side by side
  labs(x = "Counts of edits within the genic region",
       y = "Genic region",
       fill = "Infection Stage") +
  scale_fill_manual(values = c("darkmagenta", "red", "blue")) +  # Customizing colors
  theme_classic() +
  theme(legend.position = "bottom")+
  theme(axis.text.x = element_text(size = 12),  
        axis.text.y = element_text(size = 12),  
        axis.title.x = element_text(size = 16),  
        axis.title.y = element_text(size = 16),  
        legend.text = element_text(size = 14)) 

genic_region <- genic_region + plot_annotation(title = "C")

ggsave("genic_region_long.png", width = 8, height = 12)

#================== reactome pathway ============================================================================================

#  pathways with p value =< 0.05

Pre <- read.csv("Pre_pathways.csv", header = TRUE)
Pre$infection_stage <- 'Pre-infection'
Pre <- Pre %>%
  filter(Entities.pValue <= 0.05)

Mid<- read.csv("Mid_pathways.csv", header = TRUE)
Mid$infection_stage <- 'Mid-infection'
Mid <- Mid %>%
  filter(Entities.pValue <= 0.05)

Post <- read.csv("Post_pathways.csv", header = TRUE)
Post$infection_stage <- 'Post-infection'
Post <- Post %>%
  filter(Entities.pValue <= 0.05)

combined_df <- rbind(Pre, Mid, Post)
view(combined_df)

# Summarize data by Pathway and Stage
summarized_data <- combined_df %>%
  group_by(Pathway.name, infection_stage) %>%
  summarize(TotalEntities = sum(X.Entities.found))
View(summarized_data)

# Reorder the levels 
summarized_data$infection_stage <- factor(summarized_data$infection_stage, levels = c("Pre-infection", "Mid-infection", "Post-infection"))


pathway <- ggplot(summarized_data, aes(x = TotalEntities, y = Pathway.name, fill = infection_stage)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_vline(xintercept = 2, linetype = "dashed", color = "black",, size = 0.5) +  
  labs(x = "Number of genes",
       y = "Biological pathway (P value =< 0.05)",
       fill = "Infection Stage") +
  theme_classic() +
  scale_fill_manual(values = c("darkmagenta", "red", "blue")) + 
  theme_classic() +
  theme(legend.position = "bottom") +
  theme(axis.text.x = element_text(size = 12),  
        axis.text.y = element_text(size = 12),  
        axis.title.x = element_text(size = 14),  
        axis.title.y = element_text(size = 14),  
        legend.text = element_text(size = 14)) 

pathway <- pathway + plot_annotation(title = "D")

ggsave("pathway.png", width = 12, height = 10)

#===============================================================================
layout <- "
AAAAAA
BBCCDD
"
combined_plot <- Total_edits_patientwise + unique_shared + genic_region + pathway +
  plot_layout(design = layout)
