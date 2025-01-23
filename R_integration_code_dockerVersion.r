#INSTALL PACKAGES
print("Installing packages and loading libraries")
remotes::install_github("wilkelab/ggtext", upgrade = "never")

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
library(CMplot)
library(qqman)
library(ggtext)

args <- commandArgs(trailingOnly = TRUE)
keyFileDir <- "/Integratomics"

if (length(args) > 2) {
  print("Correct directory format: /path/to/input /path/to/output")
  stop("Error: too many arguments. Please check your arguments and specify input and output directories. If only one argument is provided, it will be treated as both input and output.")
} else if (length(args) == 2) {
  input_dir <- args[1]
  output_dir <- args[2]
} else if (length(args) == 1) {
  input_dir <- args[1]
  output_dir <- args[1]
} else if (length(args) == 0) {
  print("Correct directory format: /path/to/input /path/to/output")
  stop("Error: no arguments provided. Please specify input and output directories. If only one argument is provided, it will be treated as both input and output.")
}

print(paste("Input directory set to:", input_dir))
print(paste("Output directory set to:", output_dir))


#'[Settings]
workDir <- input_dir
setwd(workDir)


#'[Reading/preparing base data]
print("Reading in your data")

studyInfoDF <- tibble::tibble()
txtFiles <- list.files(paste0(workDir), pattern = "\\.txt$")
# Iterate through each .txt file and read the first 3 lines and the source file into a tibble
for (currentFile in txtFiles) {
  data <- readr::read_tsv(currentFile, n_max = 3, col_names = FALSE)
  data_t <- as_tibble(t(data)) %>%
    dplyr::filter(if_any(everything(), ~ !is.na(.) & . != ""))
  data_t <- data_t %>% mutate(source_file = currentFile)
  studyInfoDF <- bind_rows(studyInfoDF, data_t)
}
colnames(studyInfoDF) <- c("studyType", "studyName", "studyDataType", "studyFile")
#Clean up input - remove whitespaces and convert to lower
studyInfoDF <- studyInfoDF %>%
  mutate(
    studyType = tolower(gsub("\\s+", "", studyType)),
    studyDataType = tolower(gsub("\\s+", "", studyDataType))
  )
  
#Read the base file for chromosome lengths (data from UCSC)
hg38BaseFile <- read.table(paste0(keyFileDir,"/hg38_UCSC_chrom_lengths.txt"), sep="\t")

print("Preparing genome location backbone")
#Read in or prepare location backbone
#WARNING: initial preparation could take several hours because it's not paralellized and highly inefficient
if (!file.exists((paste0(keyFileDir,"/locationBackbone.txt"))))
{
  print("Preparing location backbone from hg38 base file, this could take several hours.")
  locationBackbone <- tibble::tibble()
  colnames(locationBackbone) <- c("intervalNumber","intervalChrom", "intervalStart", "intervalEnd")
  
  intervals = list()
  
  intervalCounter = 0
  step = 10000
  overlap = 5000
  
  for (i in 1:nrow(hg38BaseFile))
  {
    print(hg38BaseFile[i,])
    
    chrom <- hg38BaseFile[i,1]
    chromLength <- hg38BaseFile[i,2]
    
    for (j in seq(1, chromLength, by = step-overlap))
    {
      intervalCounter <- intervalCounter + 1
      intervalChrom <- chrom
      intervalStart <- j
      intervalEnd <- j + step -1
      
      intervals = append(intervals, c(intervalCounter, intervalChrom, intervalStart, intervalEnd))
       
    }
    
  }
  intervalsMatrix <- matrix(unlist(intervals), ncol = 4, byrow = TRUE)
  
  # Convert the matrix into a dataframe
  intervalsDF <- as.data.frame(intervalsMatrix)
  
  # Add column names
  colnames(intervalsDF) <- c("intervalNumber","intervalChrom", "intervalStart", "intervalEnd")
  
  # Append the new dataframe to the existing one
  locationBackbone <- rbind(locationBackbone, intervalsDF)
  
  #Convert to numerics
  i <- c(1, 3, 4)  
  locationBackbone[, i] <- apply(locationBackbone[, i], 2, function(x) as.numeric(as.character(x)))
  sapply(locationBackbone, class)
  
  write.table(locationBackbone, file = paste0(keyFileDir,"/locationBackbone.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
} else {
  locationBackbone <- read.table(paste0(keyFileDir,"/locationBackbone.txt"), sep="\t", header = TRUE)
}


#'[STUDY CONVERSION ------ INCOMPLETE]
#Abandoned implementation
#the idea was to accept inputs other than coordinates by automating conversion from gene symbols, probes, etc. to genomic locations


#'[TRADITIONAL (SLOW) implementation]
# for (StudyType in unique(studyInfoDF$studyType))
# {
#   fileNames <- studyInfoDF %>%
#     dplyr::filter(studyType == StudyType) %>%
#     dplyr::select(studyFile) %>%
#     as.list()
#   
#   print(StudyType)
#   print(fileNames)
#   
#   #Store all signals from the same study type into one dataframe
#   fileDataHolder <- tibble::tibble()
#   for (file in unlist(fileNames))
#   {
#     currentFileData <- readr::read_tsv(file, skip = 3, col_names = FALSE)
#     print(file)
#     fileDataHolder <- bind_rows(fileDataHolder, currentFileData)
#   }
#   colnames(fileDataHolder) <- c("intervalChrom","intervalStart", "intervalEnd", "intervalSignal")
#   fileDataHolder <- fileDataHolder %>%
#     stats::na.omit()
#   
#   #Add new column for each study type iteration
#   studyTypeColumn <- paste0("intervalSignalSum_", StudyType)
#   locationBackbone <- locationBackbone %>% mutate(!!studyTypeColumn := 0)
#   
#   # Check and sum interval signals
#   for (i in 1:nrow(fileDataHolder)) {
#     chrom <- fileDataHolder$intervalChrom[i]
#     start <- fileDataHolder$intervalStart[i]
#     end <- fileDataHolder$intervalEnd[i]
#     signal <- fileDataHolder$intervalSignal[i]
#     
#     locationBackbone <- locationBackbone %>%
#       mutate(!!studyTypeColumn := if_else(
#         intervalChrom == chrom & 
#           ((intervalStart <= start & intervalEnd >= start) | 
#              (intervalStart <= end & intervalEnd >= end)),
#         !!sym(studyTypeColumn) + signal,
#         !!sym(studyTypeColumn)
#       ))
#   }
#   
# }

#'[SUM SIGNALS, RANK SUMMED SIGNALS]
#'[OPTIMIZED]
# Optimized code
# Uses the data.table package, which works off lists
# Much faster than the traditional implementation
print("Summing signals within study types")

#Set to data.table
if (!is.data.table(locationBackbone)) {
  setDT(locationBackbone)
}

cumulativeSourceFiles <- data.table(intervalChrom = character(),
                                    intervalStart = integer(),
                                    intervalEnd = integer(),
                                    sourceFiles = character(),
                                    fileCount = integer())

for (StudyType in unique(studyInfoDF$studyType)) {
  fileNames <- studyInfoDF %>%
    dplyr::filter(studyType == StudyType) %>%
    dplyr::select(studyFile) %>%
    unlist()
  
  print(StudyType)
  print(fileNames)
  
  # Store all signals from the same study type into one dataframe
  fileDataHolder <- rbindlist(lapply(fileNames, function(file) {
    currentFileData <- fread(file, skip = 3, header = FALSE)
    setnames(currentFileData, c("intervalChrom", "intervalStart", "intervalEnd", "intervalSignal"))

    # Add the file name as a new column
    currentFileData[, sourceFile := file]

    # Ensure that currentFileData is a data.table
    if (!is.data.table(currentFileData)) {
      setDT(currentFileData)
    }

    # Convert intervalSignal to numeric
    currentFileData[, intervalSignal := ifelse(as.double(intervalSignal) == 0, .Machine$double.xmin, as.double(intervalSignal))]
    # Remove any rows where intervalStart > intervalEnd
    currentFileData <- currentFileData[intervalStart <= intervalEnd]
    
    # Transform intervalSignal to -log10(intervalSignal)
    currentFileData[, intervalSignal := -log10(intervalSignal)]

    # Set keys for both data.tables
    setkey(currentFileData, intervalChrom, intervalStart, intervalEnd)
    setkey(locationBackbone, intervalChrom, intervalStart, intervalEnd)

    # Find overlaps between currentFileData and locationBackbone
    overlaps <<- foverlaps(currentFileData, locationBackbone, type = "any", nomatch = 0)

    
    # For each interval on locationBackbone, keep only the row from currentFileData
    # with the lowest intervalSignal (p-value)
    filteredOverlaps <<- overlaps[, .SD[which.max(intervalSignal)], 
                                 by = .(intervalChrom, intervalStart, intervalEnd, sourceFile)]
    
    return(filteredOverlaps)
  }))
  
  # Add a new column for the current StudyType to store the summed interval signals
  studyTypeColumn <- paste0("intervalSignalSum_", StudyType)
  locationBackbone[, (studyTypeColumn) := 0]
  
  # Ensure fileDataHolder is a data.table
  if (!is.data.table(fileDataHolder)) {
    setDT(fileDataHolder)
  }
  
  # Sum the signals across multiple studies and update locationBackbone
  summedSignals <- fileDataHolder[, .(intervalSignalSum = sum(intervalSignal)), 
                                  by = .(intervalChrom, intervalStart, intervalEnd)]
  
  # Concatenate source files for each interval
  sourceFiles <- fileDataHolder[, .(sourceFiles = paste(unique(sourceFile), collapse = ";"),
                                    fileCount = uniqueN(sourceFile)), 
                                by = .(intervalChrom, intervalStart, intervalEnd)]
  # Update cumulativeSourceFiles with new sourceFiles information
  cumulativeSourceFiles <- rbind(cumulativeSourceFiles, sourceFiles, fill = TRUE)
  
  
  locationBackbone[sourceFiles, on = .(intervalChrom, intervalStart, intervalEnd), 
                   `:=`(sourceFiles = sourceFiles, fileCount = fileCount)]
  
  # Update the locationBackbone with the summed signals
  locationBackbone[summedSignals, on = .(intervalChrom, intervalStart, intervalEnd), 
                   (studyTypeColumn) := intervalSignalSum]
  
  # Add a ranking column for the current StudyType
  rankColumn <- paste0("rank_", StudyType)
  locationBackbone[, (rankColumn) := rank(-get(studyTypeColumn), ties.method = "max")]
  
  # Add a column to store the concatenated source file information
  locationBackbone[sourceFiles, on = .(intervalChrom, intervalStart, intervalEnd), 
                   sourceFiles := sourceFiles]
}

cumulativeSourceFiles <- cumulativeSourceFiles[, .(sourceFiles = paste(unique(sourceFiles), collapse = ";"),
                                                   fileCount = sum(fileCount)), 
                                               by = .(intervalChrom, intervalStart, intervalEnd)]
#Add source file data to the backbone
locationBackbone[cumulativeSourceFiles, on = .(intervalChrom, intervalStart, intervalEnd), 
                 `:=`(sourceFiles = sourceFiles, fileCount = fileCount)]

#Convert to tibbles, update locationBackbone, return to data.table
locationBackbone_tbl <- as_tibble(locationBackbone)
cumulativeSourceFiles_tbl <- as_tibble(cumulativeSourceFiles)


locationBackbone_tbl <- locationBackbone_tbl %>%
  dplyr::left_join(cumulativeSourceFiles_tbl, by = c("intervalChrom", "intervalStart", "intervalEnd")) %>%
  dplyr::mutate(sourceFiles = coalesce(sourceFiles.y, sourceFiles.x),
         fileCount = coalesce(fileCount.y, fileCount.x)) %>%
  dplyr::select(-sourceFiles.x, -sourceFiles.y, -fileCount.x, -fileCount.y)

locationBackbone <- as.data.table(locationBackbone_tbl) %>%
  dplyr::mutate_at(c('fileCount'), ~replace_na(.,0)) %>%
  dplyr::mutate(sourceFiles = ifelse(is.na(sourceFiles), "", sourceFiles))

#Add information on how many signals are actually mapped from each individual file - helps to double-check your data
print("Checking number of signals mapped to intervals from each individual file")
uniqueSourceFiles <- unique(unlist(strsplit(locationBackbone_tbl$sourceFiles, ";"))) %>%
  as.data.frame(nm = (c("sourceFiles"))) %>%
  dplyr::filter(!is.na(sourceFiles))

uniqueSourceFiles$numIntervals <- sapply(uniqueSourceFiles$sourceFiles, function(sourceFile) {
  sum(str_count(locationBackbone$sourceFiles, fixed(sourceFile)))
})

#'[ADD ARITHMETIC MEAN OF RANKS - rank product]
#You can use either geometric or arithmetic mean (Breitling et al. 2016), we used the arithmetic mean
# Identify columns whose name contains "rank_", apply arithmetic mean across intervals
# arithm_rank_product is therefore the rank product
study_rank_cols <- grep("rank_", names(locationBackbone), value = TRUE)
#geo_mean <- function(x) {
#  exp(mean(log(x)))
#}
locationBackbone[, interval_rank_product := apply(.SD, 1, mean), .SDcols = study_rank_cols]


# Add a ranking column where ties are assigned the highest rank
# Example : if three intervals, A, B and C, are tied for 1st place, they will each be assigned the rank "3"
locationBackbone[, arithm_rank_product := frank(interval_rank_product, ties.method = "max")]

#'[GENE DENSITY ESTIMATION]
print("Estimating gene densities within intervals")
#Find the locations of intervals with the same gene density, for the data on gene density we will use the Ensembl GRCh38.p14 assembly
#We included protein-coding genes from chr 1-22+X+Y, that have a UCSC Gene Stable ID (20 037 genes)

#Explanation: permuting within intervals that have the same gene density prevents bias based on how gene-rich a certain part of the genome is
#For example, if you permuted signals from intergenic regions with gene regions, it would significantly reduce the threshold for a signal being statistically significant
#Thus, permuting intervals together based on their gene density is done to avoid flooding the results with false positives

biomartGeneLocations <- read.table(paste0(keyFileDir,"/mart_export.txt"), sep="\t", header = TRUE) %>%
  dplyr::distinct() %>%
  dplyr::select(Chromosome.scaffold.name, Gene.start..bp., Gene.end..bp., Gene.name, Gene.stable.ID) %>%
  dplyr::rename(Chrom = Chromosome.scaffold.name, Start = Gene.start..bp., End = Gene.end..bp., Gene_name = Gene.name, Ensembl_ID = Gene.stable.ID) %>%
  data.table::as.data.table()

biomartGeneLocations[, Chrom := paste0("chr", Chrom)]


# Convert biomartGeneLocations to a data.table
setDT(biomartGeneLocations)

# Set keys for interval join
setkey(locationBackbone, intervalChrom, intervalStart, intervalEnd)
setkey(biomartGeneLocations, Chrom, Start, End)

# Perform the interval join and count the number of genes in each interval
overlaps <- foverlaps(biomartGeneLocations,
                      locationBackbone,
                      by.x = c("Chrom",
                               "Start",
                               "End"),
                      by.y = c("intervalChrom",
                               "intervalStart",
                               "intervalEnd"),
                      type = "any",
                      nomatch = 0)
colnames(overlaps)[colnames(overlaps) == "Chrom"] = "intervalChrom"

# Count the number of genes and concatenate gene names for each interval
geneCounts <- overlaps[, .(
  gene_count = .N,
  genes = paste(unique(Gene_name), collapse = ";"),
  geneEnsemblIDs = paste(unique(Ensembl_ID), collapse = ";")
), by = .(intervalChrom, intervalStart, intervalEnd)]

# Add gene count and concatenated gene names to locationBackbone
locationBackbone <- locationBackbone[
  geneCounts,
  on = .(intervalChrom, intervalStart, intervalEnd),
  `:=`(
    gene_count = i.gene_count,
    gene_names = i.genes,
    gene_Ensembl_IDs = i.geneEnsemblIDs
  )
]

# Ensure that intervals with no genes have NA or zero
locationBackbone[is.na(gene_count), `:=`(gene_count = 0, gene_names = "")]
#Remove the intermediates "genes" column
locationBackbone[, genes := NULL]
#Cleanup
rm(overlaps)
rm(geneCounts)

#Number of intervals with a specific number of genes
#print(table(locationBackbone$gene_count))

#'[REMOVE INTERVALS WITH NO SIGNALS]

print("Removing intervals with no signals")
#For the analysis we only need locations that have any signals
#This is done at this stage so that merging singleton gene density regions can be done accurately (see next section)

#Remove unnecessary columns, keep only intervals with signal (if arithm_rank_product is the same as the lowest rank (length of data), then it's signal-less)
#Data cleanup - Convert all rank_ columns to numerics
nonZeroLocations <- locationBackbone %>%
  dplyr::select(intervalNumber, intervalChrom, intervalStart, intervalEnd, sourceFiles, fileCount, interval_rank_product, arithm_rank_product, gene_count, gene_names, gene_Ensembl_IDs, dplyr::starts_with("rank_")) %>%
  dplyr::filter(arithm_rank_product != length(locationBackbone$intervalNumber)) %>%
  dplyr::mutate(across(starts_with("rank_"), as.numeric))

#'[MERGE GENE DENSITIES]
print("Merging gene densities for accurate permutations")
#### NOTE ON GENE DENSITIES
# Certain regions have very high gene densities because they are populated by a multitude of small genes
# An example is an interval on chr5 (141510001 - 141520000), which has 23 genes on it, and is the only one to have that many
# This poses a problem for permutations, as intervals with a specific gene count will get permuted with other intervals with the same gene count
# This is done in order to account for signal bias based on gene density
# However, you will note that if you have only 1 interval with a given gene density, it will be permuted ... with itself
# This could cause wildly inaccurate results later on down the line in specific outlier cases
# Therefore, we will group singletons together until they can form groups big enough to feasibly be permuted at least 1000 times with different results
# (1000 permutations is the default value for the number of permutations)
minimumDensityGroupSize <- 1000

# Count the number different study types - columns that start with "rank_"
rank_columns_count <- sum(grepl("^rank_", colnames(nonZeroLocations)))
# Get unique gene_count values in descending order
unique_gene_counts <- sort(unique(nonZeroLocations$gene_count), decreasing = TRUE)
# Calculate the number of combinations (A!^B), where A is the number of intervals with the same gene count and B is number of study types
# If the number of combinations isn't sufficient, merge down to the next gene_count group
for (current_gene_count in unique_gene_counts) {
  # Calculate the factorial of the count of intervals with this gene_count
  count_intervals <- sum(nonZeroLocations$gene_count == current_gene_count)
  factorial_value <- factorial(count_intervals)
  threshold_permutation <- factorial_value^rank_columns_count
  # Check if the factorial value is less than minimumDensityGroupSize
  if (threshold_permutation < minimumDensityGroupSize) {
    # Decrement the gene_count for these intervals
    nonZeroLocations[gene_count == current_gene_count, gene_count := gene_count - 1]
  }
}
unique_gene_counts_merged <- sort(unique(nonZeroLocations$gene_count), decreasing = TRUE)
print("The following gene density regions were kept after merging:")
print(unique_gene_counts_merged)
#print(table(nonZeroLocations$gene_count))

#'[PERMUTATIONS]

#We're using a method based on Breitling et al. 2004 (https://doi.org/10.1016/j.febslet.2004.07.055)
#This method is more rigorous, but also more computationally intensive
#The are 2 minor differences between our method and Breitling's: 
#1. summing of signals within a single data source (study type)
#2. permutations take place within intervals of the same gene density to account for disparity between regions

#Note: You could estimate p-values by fitting a distribution to your data
#We've tried this, but no distribution fits particularly well (proven by goodness of fit tests), leading to inaccurate estimates

###############################
### OPTIMIZATION VERSION ### V12 - TIDYR
###############################
nPerm <- 1000



tibble_rankSum <- function(ranks, numPerm) {
  tibble::as_tibble(ranks, .name_repair = "unique_quiet") %>% 
    dplyr::mutate(
      RS = rowMeans(dplyr::across(tidyselect::starts_with("rank_")), na.rm=TRUE)
    ) %>% 
    tidyr::pivot_longer(tidyselect::starts_with("rank_"), names_to = "study", values_to = "rank") %>% 
    dplyr::mutate(across(starts_with("rank_"), as.integer)) %>% 
    dplyr::arrange(study) %>% 
    dplyr::bind_cols(matrix(data=NA_integer_, nrow=nrow(.), ncol=numPerm, dimnames=list(NULL, paste0("RSperm",1:numPerm)))) %>% 
    tidyr::pivot_longer(tidyselect::matches("RSperm"), names_to = "perm", values_to = "RSperm") %>% 
    dplyr::select(-RSperm) %>% 
    dplyr::mutate(perm = stringr::str_remove(perm, "RSperm") %>% as.integer()) %>% 
    ## perms
    dplyr::group_by(study, perm) %>% 
    dplyr::mutate(sampleRank = sample(rank)) %>% 
    dplyr::group_by(intervalNumber, perm) %>% 
    dplyr::left_join(dplyr::reframe(., RSperm = mean(sampleRank)), by=c("intervalNumber", "perm")) %>%
    ## summarise
    tidyr::pivot_wider(id_cols = c(intervalNumber, study, RS), names_from = perm, names_prefix="RSperm", values_from = RSperm) %>% 
    dplyr::summarise(dplyr::across(tidyselect::starts_with("RS"), unique))
}
start <- Sys.time()
print("!----------------------------------------------------------------------------!")
print("Beginning permutation step - this is the most computationally intensive step")
print("For ~160k intervals, 5 groups and 1000 permutations the expected time is ~1.5 hours and peaks at ~160GB of RAM")
print("WARNING: MONITOR YOUR RAM USAGE, THIS STEP CAN BE VERY MEMORY INTENSIVE")
print(paste0("Current setting: ", nPerm, " permutations."))
print(paste0("Started at: ", start))
print("!----------------------------------------------------------------------------!")

result <- nonZeroLocations %>%
  dplyr::group_by(gene_count) %>%
  dplyr::group_modify(~ tibble_rankSum(.x, numPerm = nPerm)) %>%
  dplyr::ungroup()

print(paste0("Permutations completed. Time needed: ", Sys.time() - start))

#'[STATISTICS]
#Add row that counts how many values in the permuted RP tibble are <= the real RP value (RPperm <= RPreal)
#We will call this the cValue
print("Running statistics on permuted ranks - this also takes some time")
print(paste0("Process started at: ", Sys.time()))
start <- Sys.time()

print("Calculating c-value")
result2 <- result %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    cValue = sum(RS >= c_across(starts_with("RSPerm")))
  ) %>%
  dplyr::ungroup()
print(paste0("Process completed. Time: ", Sys.time() - start))
print("Calculating empirical p-values")

#Adding pseudocount
result2 <- result2 %>%
  dplyr::mutate(cValuePseudoCount = cValue + 1)
#Store backup
resultRaw <- result2

#Append location data
locationData <- nonZeroLocations %>%
  dplyr::select(intervalNumber, intervalChrom, intervalStart, intervalEnd, sourceFiles, fileCount, gene_count, gene_names, gene_Ensembl_IDs, arithm_rank_product)
result2 <- result2 %>%
  dplyr::select(!tidyselect::starts_with("RSPerm"))
result2 <- dplyr::left_join(result2, locationData, by = "intervalNumber") %>%
  dplyr::mutate(arithm_rank_desc = rank(-RS, ties.method = "max"))

#Add row that calculates the expected RP value (Erp = cValue/number of permutations)
#Add row that calculates the p-value, in this case the percentage of false positives (PFP), where PFP = Erp/rank

result2 <- result2 %>%
  dplyr::mutate(empiricalPval = (cValuePseudoCount) / (nPerm+1)) %>%
  dplyr::mutate(empiricalPval = round(empiricalPval, digits = 3))



print("Statistics calculated, outputting result files")

#Output raw results
results_for_print_raw <- result2 %>%
  dplyr::rename(rank_product = RS) %>%
  dplyr::rename(interval_rank = arithm_rank_desc) %>%
  dplyr::select(intervalNumber, intervalStart, intervalEnd, sourceFiles, fileCount, gene_count, gene_names, gene_Ensembl_IDs, interval_rank, rank_product, cValue, cValuePseudoCount, empiricalPval)
write.table(results_for_print_raw, file=paste0(output_dir, "Integration_results_preBoot_raw_", Sys.Date(), ".tsv"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

#Output significant results
results_for_print_sig <- results_for_print_raw %>%
  dplyr::filter(empiricalPval <= 0.05)
write.table(results_for_print_sig, file=paste0(output_dir, "Integration_results_preBoot_significant_", Sys.Date(), ".tsv"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

#Output unique gene names mapped to significant intervals
sigGenes <- results_for_print_sig %>%
  dplyr::select(gene_names) %>%
  tidyr::separate_longer_delim(gene_names, delim = ";") %>%
  dplyr::distinct() %>%
  dplyr::filter(gene_names != "", !is.na(gene_names))
write.table(sigGenes, file=paste0(output_dir, "Genes_from_preBoot_significant_intervals_", Sys.Date(), ".tsv"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

sigEnsembl <- results_for_print_sig %>%
  dplyr::select(gene_Ensembl_IDs) %>%
  tidyr::separate_longer_delim(gene_Ensembl_IDs, delim = ";") %>%
  dplyr::distinct() %>%
  dplyr::filter(gene_Ensembl_IDs != "", !is.na(gene_Ensembl_IDs))
write.table(sigEnsembl, file=paste0(output_dir, "EnsemblIDs_from_preBoot_significant_intervals_", Sys.Date(), ".tsv"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

write.table(uniqueSourceFiles, file=paste0(output_dir, "Signals_by_file_", Sys.Date(), ".tsv"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

#'[PLOT OUTPUTS]
print("Creating P-value histogram and various Manhattan plots")

png(filename = paste0(output_dir, "empiricalPvalueHistogram.png"),
    width = 1280,
    height = 840)
ggplot(data = result2, aes(x=empiricalPval)) + 
  geom_histogram(binwidth=1/nPerm)
dev.off()

#_____________________
#Making Manhattan plots
#_____________________
#Add -log10(pval)
results_for_print_raw$NegLog10Signal <- -log10(results_for_print_raw$empiricalPval)

#Sort chromosomes numerically and assign positions for x-axis, define threshold for coloring columns red above the Y-axis line
results_for_print_raw$intervalChrom <- factor(results_for_print_raw$intervalChrom, levels = unique(results_for_print_raw$intervalChrom))
results_for_print_raw$Position <- as.numeric(results_for_print_raw$intervalChrom) + (results_for_print_raw$intervalStart / max(results_for_print_raw$intervalStart))
threshold <- -log10(0.05)

chromosome_midpoints <- results_for_print_raw %>%
  group_by(intervalChrom) %>%
  summarize(mid_point = mean(Position))

png(filename = paste0(output_dir, "Integ_Manhattan_plot_preBoot_v1_", Sys.Date(), ".png"),
    width = 1600,
    height = 840)
ggplot(results_for_print_raw, aes(x = Position, y = NegLog10Signal)) +
  geom_point(aes(color = NegLog10Signal > threshold), alpha = 0.7) +
  scale_color_manual(values = c("black", "red"), guide = "none") +
  geom_hline(yintercept = threshold, linetype = "dashed", color = "blue") +
  labs(
    title = "Manhattan Plot of Genomic Regions",
    x = "Chromosomes",
    y = expression(-log[10](empirical_pValue))
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

Manhattan_results <- results_for_print_raw %>%
  dplyr::mutate(intervalChrom = recode(intervalChrom,
                             "chrX" = "chr23",
                             "chrY" = "chr24",
                             "chrM" = "chr25")) %>%
  dplyr::mutate(across('intervalChrom', str_replace, 'chr', '')) %>%
  dplyr::mutate(intervalChrom = as.numeric(intervalChrom))

chromosome_midpoints2 <- Manhattan_results %>%
  group_by(intervalChrom) %>%
  summarize(mid_point = mean(Position))

png(filename = paste0(output_dir, "Integ_Manhattan_plot_preBoot_v2_", Sys.Date(), ".png"),
    width = 1600,
    height = 840) 
ggplot(Manhattan_results, aes(x = Position, y = NegLog10Signal)) +
  geom_point( aes(color=as.factor(intervalChrom)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
  geom_hline(yintercept = threshold, linetype = "dashed", color = "red") +
  scale_x_continuous(
    label = chromosome_midpoints2$intervalChrom,
    breaks= chromosome_midpoints2$mid_point
    ) +
  scale_y_continuous(expand = c(0, 0) ) +
  
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )
dev.off()

png(filename = paste0(output_dir, "Integ_Manhattan_plot_preBoot_v3_", Sys.Date(), ".png"),
    width = 1600,
    height = 840) 
qqman::manhattan(Manhattan_results, chr="intervalChrom", bp="intervalStart", snp="intervalNumber", p="empiricalPval", genomewideline = -log10(0.05) )
dev.off()


sigIntervals <- Manhattan_results %>%
  dplyr::filter(empiricalPval <= 0.05) %>%
  dplyr::select(intervalNumber) %>%
  unlist() %>%
  as.vector()
png(filename = paste0(output_dir, "Integ_Manhattan_plot_preBoot_v4_", Sys.Date(), ".png"),
    width = 1600,
    height = 840) 
qqman::manhattan(Manhattan_results, chr="intervalChrom", bp="intervalStart", snp="intervalNumber", p="empiricalPval", genomewideline = -log10(0.05), highlight = sigIntervals)
dev.off()


Manhattan_results_circ <- Manhattan_results %>%
  dplyr::select(intervalNumber, intervalChrom, intervalStart, empiricalPval)
  
CMplot(Manhattan_results_circ,
       plot.type="m",
       col=c("grey30","grey60"),
       LOG10=TRUE,
       threshold=c(0.05),
       threshold.lty=c(1,2),
       threshold.lwd=c(1,1),
       threshold.col=c("black","grey"),
       amplify=TRUE,
       chr.den.col=NULL,
       signal.col=c("red"),
       signal.cex=1,
       signal.pch=c(19,19),
       file.name = paste0(output_dir, "Integ_Manhattan_plot_preBoot_v5_", Sys.Date()),
       file="jpg",
       file.output=TRUE,
       verbose=TRUE,
       dpi=300,
       width=14,height=6)

CMplot(Manhattan_results_circ, plot.type="c", r=1.6,
       outward=TRUE, cir.chr.h=.1 ,chr.den.col="orange",
       file.name = paste0(output_dir, "Integ_Manhattan_plot_preBoot_v6_", Sys.Date()),
       file="jpg",
       dpi=300, chr.labels=seq(1,25))

CMplot(
  Manhattan_results_circ,
  plot.type = "c",
  r = 1.6,
  outward = FALSE,
  col=c("grey30","grey60"),
  cir.chr.h = 1,
  signal.cex = 0.1,
  file.name = paste0(output_dir, "Integ_Manhattan_plot_preBoot_v7_", Sys.Date()),
  file = "jpg",
  dpi = 300,
  chr.labels = seq(1, 25)
)

CMplot(
  Manhattan_results_circ,
  plot.type = "c",
  r = 1.6,
  outward = TRUE,
  col=c("grey30","grey60"),
  cir.chr.h = 1,
  signal.cex = 0.1,
  file.name = paste0(output_dir, "Integ_Manhattan_plot_preBoot_v8_", Sys.Date()),
  file = "jpg",
  dpi = 300,
  chr.labels = seq(1, 25)
)


print("Step 1 of analysis complete!")

if (format(Sys.Date(), "%m") == "01" && as.integer(format(Sys.Date(), "%d")) >= 1 && as.integer(format(Sys.Date(), "%d")) <= 5) {
  print("And have a happy new year!")
}

