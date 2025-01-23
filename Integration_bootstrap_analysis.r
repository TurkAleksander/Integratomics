#INSTALL PACKAGES
print("Installing packages and loading libraries")

install.packages("dplyr")
install.packages("ggplot2")
install.packages("readr")
install.packages("stringr")
install.packages("tibble")
install.packages("data.table")
install.packages('MASS')
install.packages("tidyr")
install.packages("purrr")
install.packages("furrr")
install.packages("future")
install.packages("ggtext")
install.packages("qqman")
install.packages("CMplot")

library(CMplot)
library(qqman)
library(dplyr)
library(ggplot2)
library(readr)
library(stringr)
library(tibble)
library(data.table)
library(MASS)
library(tidyr)
library(purrr)
library(furrr)
library(future)
library(ggtext)



#'[Settings]
workDir <- "/cmg1scratch/PROJECTS/MS_integratomics/Data_integration/All_integ_data"
setwd(workDir)

#'[Reading in the data]

print("Reading in real dataset result file")
fileTarget <- list.files(workDir, pattern = "^Integration_results_raw_real", full.names = TRUE, ignore.case = TRUE)
realData <-  read.table(fileTarget, header = TRUE, sep = "\t")

#Find the bootstrapped datasets
bootstrappedFiles <- list.files(workDir, pattern = "^Integration_results_raw_bootstrapRun", full.names = TRUE)
key_column <- "intervalNumber"

print("Appending data from bootstrapped result files.")
print("This takes a few minutes ...")
mergedData <- realData
#Loop through each bootstrapped file
for (i in seq_along(bootstrappedFiles)) {
  #Read the file, check for column
  boot_data <- read.table(bootstrappedFiles[i], header = TRUE, sep = "\t")
  if (!all(c(key_column, "empiricalPval") %in% colnames(boot_data))) {
    stop(paste("File", bootstrappedFiles[i], "does not contain the required columns."))
  }
  
  #Rename the empiricalPval column for each file
  col_name <- paste0("boot_empiricalPval_", i)
  boot_data <- boot_data %>%
    dplyr::select(all_of(key_column), empiricalPval) %>%
    dplyr::rename(!!col_name := empiricalPval)
  
  # Join the current bootstrapped data with realData
  mergedData <- mergedData %>%
    dplyr::left_join(boot_data, by = key_column)
}

#'[Boot statistics from distribution]
#We are essentially checking the distribution of the bootstrapped p-values
#We're doing this by checking each genomic interval's confidence interval through the mean and standard deviation (SD)
#If the 95% confidence interval (CI) (mean + 1.96*SD) is lower than our alpha (0.05), then this is a good genomic interval

#Old slower implementation
# mergedData2 <- mergedData %>%
#   dplyr::rowwise() %>%
#   dplyr::mutate(
#     boot_avg = mean(c_across(starts_with("boot_empiricalPval")), na.rm =TRUE),
#     boot_stdev = sd(c_across(starts_with("boot_empiricalPval")),  na.rm =TRUE)) %>%
#   dplyr::ungroup()


print("Calculating bootstrapped P-value")
print("This will be done quickly ...")

mergedData2 <- mergedData %>%
  dplyr::mutate(
    boot_avg = rowMeans(dplyr::select(., starts_with("boot_empiricalPval")), na.rm = TRUE),
    boot_stdev = apply(dplyr::select(., starts_with("boot_empiricalPval")), 1, sd, na.rm = TRUE)) %>%
  dplyr::mutate(bootstrap_dist = (boot_avg + 1.96*boot_stdev))
  
print("Calculated. Outputting results file")



##'[Boot statistics [old]]
# Assuming mergedData is your dataframe
# mergedData <- mergedData %>%
#   dplyr::rowwise() %>%
#   dplyr::mutate(
#     count_lower = sum(c_across(starts_with("boot_empiricalPval")) < empiricalPval, na.rm = TRUE),
#     non_empty_cells = sum(!is.na(c_across(starts_with("boot_empiricalPval")))),
#     bootstrap_reliability_Pval = (count_lower + 1) / (non_empty_cells + 1)
#   ) %>%
#   dplyr::ungroup() %>%
#   dplyr::select(-count_lower, -non_empty_cells) %>%
#   dplyr::mutate(bootstrap_reliability_Pval = round(bootstrap_reliability_Pval, digits = 3))

mergedData2 <- mergedData2 %>%
  dplyr::select(-starts_with("boot_empiricalPval"))
write.table(mergedData2, file=paste0("Integration_bootstrappedResults_raw_", Sys.Date(), ".tsv"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

print("Outputting graphs")

png(filename = "bootstrappingDistHistogram.png",
    width = 1280,
    height = 840)
ggplot(data = mergedData2, aes(x=bootstrap_dist)) + 
  geom_histogram(binwidth=1/1000)
dev.off()

#_____________________
#Making Manhattan plots
#_____________________
#Add -log10(pval)
mergedData2$NegLog10Signal <- -log10(mergedData2$bootstrap_dist)
mergedData2 <- mergedData2 %>%
  dplyr::mutate(NegLog10Signal = replace(NegLog10Signal, NegLog10Signal < 0, 0))

#Sort chromosomes numerically and assign positions for x-axis, define threshold for coloring columns red above the Y-axis line
mergedData2$intervalChrom <- factor(mergedData2$intervalChrom, levels = unique(mergedData2$intervalChrom))
mergedData2$Position <- as.numeric(mergedData2$intervalChrom) + (mergedData2$intervalStart / max(mergedData2$intervalStart))
threshold <- -log10(0.05)

chromosome_midpoints <- mergedData2 %>%
  dplyr::group_by(intervalChrom) %>%
  dplyr::summarize(mid_point = mean(Position))

png(filename = paste0("Integ_bootstrap_Manhattan_plot_v1_", Sys.Date(), ".png"),
    width = 1600,
    height = 840)
ggplot(mergedData2, aes(x = Position, y = NegLog10Signal)) +
  geom_point(aes(color = NegLog10Signal > threshold), alpha = 0.7) +
  scale_color_manual(values = c("black", "red"), guide = "none") +
  geom_hline(yintercept = threshold, linetype = "dashed", color = "blue") +
  labs(
    title = "Manhattan Plot of Genomic Regions",
    x = "Chromosomes",
    y = expression(-log[10](bootstrap_dist))
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  ) +
  scale_x_continuous(
    breaks = chromosome_midpoints$mid_point,
    labels = chromosome_midpoints$intervalChrom
  )
dev.off()

