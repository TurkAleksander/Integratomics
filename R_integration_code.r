#INSTALL PACKAGES
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
#install.packages('VGAM')


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
#library(VGAM)


#'[Settings]
workDir <- "/cmg1scratch/PROJECTS/MS_integratomics/Data_integration/TestData"
setwd(workDir)


#'[Reading/preparing base data]
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
hg38BaseFile <- read.table("/cmg1scratch/PROJECTS/MS_integratomics/Data_integration/hg38_UCSC_chrom_lengths.txt", sep="\t")


#Read in or prepare location backbone
#WARNING: initial preparation could take several hours because it's not paralellized and highly inefficient
if (!file.exists("/cmg1scratch/PROJECTS/MS_integratomics/Data_integration/hg38_UCSC_chrom_lengths.txt"))
{
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
  
  # Add column names (optional)
  colnames(intervalsDF) <- c("intervalNumber","intervalChrom", "intervalStart", "intervalEnd")
  
  # Append the new dataframe to the existing one
  locationBackbone <- rbind(locationBackbone, intervalsDF)
  
  #Convert to numerics
  i <- c(1, 3, 4)  
  locationBackbone[, i] <- apply(locationBackbone[, i], 2, function(x) as.numeric(as.character(x)))
  sapply(locationBackbone, class)
  
  write.table(locationBackbone, file = "locationBackbone_R-version.txt", sep = "\t", row.names = FALSE, quote = FALSE)
} else {
  locationBackbone <- read.table("/cmg1scratch/PROJECTS/MS_integratomics/Data_integration/locationBackbone.txt", sep="\t", header = TRUE)
}


#'[STUDY CONVERSION ------ INCOMPLETE]
#At this stage we go through the headers of studies in the directory
#The goal is to convert all data into coordinates and translate it to hg38
#By default it assumes that coordinates are hg38 if nothing is specified

for (studyIndex in 1:nrow(studyInfoDF)){
  
  #Check coordinates
  if ("coor" %in% tolower(studyInfoDF[studyIndex,3])){
    #print(studyIndex)
    
    if (("co") %in% tolower(studyInfoDF[studyIndex,3]) | ("or") %in% tolower(studyInfoDF[studyIndex,3])) {
      print(studyIndex)
    }
    
  }
}




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
#Copilot code optimization
#Uses the data.table package, which works off lists
#Much faster, but I don't fully understand the functions. Results seem to be identical, though
setDT(locationBackbone)

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
    return(currentFileData)
  }))
  
  fileDataHolder <- fileDataHolder %>%
    na.omit() %>%
    .[intervalStart < intervalEnd]
  
  # Add a new column for the current StudyType to store the summed interval signals
  studyTypeColumn <- paste0("intervalSignalSum_", StudyType)
  locationBackbone[, (studyTypeColumn) := 0]
  
  # Convert intervalSignal to numeric
  fileDataHolder[, intervalSignal := as.numeric(intervalSignal)]
  
  # Perform the interval join and sum the signals
  setkey(fileDataHolder, intervalChrom, intervalStart, intervalEnd)
  setkey(locationBackbone, intervalChrom, intervalStart, intervalEnd)
  
  overlaps <- foverlaps(fileDataHolder, locationBackbone, type = "any", nomatch = 0)
  summedSignals <- overlaps[, .(intervalSignalSum = sum(intervalSignal)), 
                            by = .(intervalChrom, intervalStart, intervalEnd)]
  
  # Update the locationBackbone with the summed signals
  locationBackbone[summedSignals, on = .(intervalChrom, intervalStart, intervalEnd), 
                   (studyTypeColumn) := intervalSignalSum]
  
  # Add a ranking column for the current StudyType
  rankColumn <- paste0("rank_", StudyType)
  locationBackbone[, (rankColumn) := rank(-get(studyTypeColumn), ties.method = "max")]
}

#'[ADD GEOMETRIC MEAN OF RANKS]
# Identify columns whose name contains "rank_", apply geometric mean across intervals
study_rank_cols <- grep("rank_", names(locationBackbone), value = TRUE)
geo_mean <- function(x) {
  exp(mean(log(x)))
}
locationBackbone[, interval_rank_product := apply(.SD, 1, geo_mean), .SDcols = study_rank_cols]


# Add a ranking column where ties are assigned the highest rank
# Example : if intervals A, B and C are tied for 1st place, they will each be assigned the rank "3"
locationBackbone[, geo_mean_rank := frank(interval_rank_product, ties.method = "max")]

#'[GENE DENSITY ESTIMATION]

#Find the locations of intervals with the same gene density, for the data on gene density we will use the Ensembl GRCh38.p14 assembly
#We included protein-coding genes from chr 1-22+X+Y+MT, that have a UCSC Gene Stable ID (20 037 genes)

#Explanation: permuting within intervals that have the same gene density prevents bias based on how gene-rich a certain part of the genome is
#For example, if you permuted signals from intergenic regions with gene regions, it would significantly reduce the threshold for a signal being statistically significant
#Thus, permuting intervals together based on their gene density is done to avoid flooding the results with false positives

biomartGeneLocations <- read.table("/cmg1scratch/PROJECTS/MS_integratomics/Data_integration/mart_export.txt", sep="\t", header = TRUE) %>%
  dplyr::distinct() %>%
  dplyr::select(Chromosome.scaffold.name, Gene.start..bp., Gene.end..bp., Gene.name, Gene.stable.ID) %>%
  dplyr::rename(Chrom = Chromosome.scaffold.name, Start = Gene.start..bp., End = Gene.end..bp., Gene_name = Gene.name, Ensembl_ID = Gene.stable.ID) %>%
  data.table::as.data.table()

biomartGeneLocations[, Chrom := paste0("chr", Chrom)]


# Ensure biomartGeneLocations is a data.table
setDT(biomartGeneLocations)

# Set keys for interval join
setkey(locationBackbone, intervalChrom, intervalStart, intervalEnd)
setkey(biomartGeneLocations, Chrom, Start, End)

# Perform the interval join and count the number of genes in each interval
overlaps <- foverlaps(biomartGeneLocations, locationBackbone, by.x = c("Chrom", "Start", "End"), by.y = c("intervalChrom", "intervalStart", "intervalEnd"), type = "any", nomatch = 0)
colnames(overlaps)[colnames(overlaps) == "Chrom"] = "intervalChrom"
geneCounts <- overlaps[, .N, by = .(intervalChrom, intervalStart, intervalEnd)]

# Add the count of genes to locationBackbone
locationBackbone[, gene_count := 0]
locationBackbone[geneCounts, on = .(intervalChrom, intervalStart, intervalEnd), gene_count := N]

#Cleanup
rm(overlaps)
rm(geneCounts)

#Number of intervals with a specific number of genes
print(table(locationBackbone$gene_count))

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
rank_columns_count <- sum(grepl("^rank_", colnames(locationBackbone)))
# Get unique gene_count values in descending order
unique_gene_counts <- sort(unique(locationBackbone$gene_count), decreasing = TRUE)
print(unique_gene_counts)
# Calculate the number of combinations (A!^B), where A is the number of intervals with the same gene count and B is number of study types
# If the number of combinations isn't sufficient, merge down to the next gene_count group
for (current_gene_count in unique_gene_counts) {
  # Calculate the factorial of the count of intervals with this gene_count
  count_intervals <- sum(locationBackbone$gene_count == current_gene_count)
  factorial_value <- factorial(count_intervals)
  threshold_permutation <- factorial_value^rank_columns_count
  # Check if the factorial value is less than minimumDensityGroupSize
  if (threshold_permutation < minimumDensityGroupSize) {
    # Decrement the gene_count for these intervals
    locationBackbone[gene_count == current_gene_count, gene_count := gene_count - 1]
  }
}

print(table(locationBackbone$gene_count))

#'[PERMUTATIONS]
#
#Note: You could estimate p-values by fitting a distribution to your data
#We've tried this, but no distribution fits particularly well (proven by goodness of fit tests), leading to inaccurate estimates
#Thus, we're using a method based on Breitling et al. 2004 (https://doi.org/10.1016/j.febslet.2004.07.055)
#This method is more rigorous, but also more computationally intensive
#The are 2 minor differences between our method and Breitling's: 
#1. summing of signals within a single data source (study type)
#2. permutations take place within intervals of the same gene density to account for disparity between regions


###############################
### TROUBLESHOOTING VERSION ### V2!!!!
###############################
# Temporarily disable parallel processing for debugging

start <- Sys.time()

if (.Platform$OS.type == "unix"){
  plan(multicore)
} else if  (.Platform$OS.type == "windows") {
  plan(multisession)
} else {
  plan(sequential)
}

calculate_permuted_rank_products <- function(data, num_permutations = 1000) {
  #Convert input data to data.table format for efficiency
  setDT(data)
  
  #Define function for calculating geometric mean - used later
  geometric_mean_rank <- function(ranks) {
    ranks <- as.numeric(ranks)  # Explicitly convert to numeric
    if (any(is.na(ranks))) {
      stop("Ranks contain NA values")
    }
    exp(mean(log(ranks)))
  }
  
  # Iterate over each intervalNumber in the data and perform permutations sequentially
  results_list <- map(data$intervalNumber, function(interval_num) {
    if (interval_num*10 %% 5000 <= 10) {
      print(paste("Processing interval:", interval_num))  # Debug print
    }
    #print(paste0("Permutation section starts - current interval: ", interval_num))
    current_interval <- data[data$intervalNumber == interval_num, ]
    same_gene_count_intervals <- data[data$gene_count == current_interval$gene_count, ]
    #print(paste0("Printing same_gene_count_intervals head:", head(same_gene_count_intervals)))
      
    permuted_rank_products <- replicate(num_permutations, {
      #message("Checkpoint 1 inside permuted_rank_products")
      
      # Convert same_gene_count_intervals to data.table
      permuted_ranks <- as.data.table(same_gene_count_intervals)
      #message("Checkpoint 2 inside permuted_rank_products")
      
      # Sample the ranks in each group
      permuted_ranks <- permuted_ranks[, lapply(.SD, function(x) sample(x, .N, replace = FALSE)), 
                                       by = intervalNumber, .SDcols = patterns("^rank_")]
      
      # Set key for easy subsetting by intervalNumber
      setkey(permuted_ranks, intervalNumber)
      #message("Checkpoint 3 inside permuted_rank_products")
      
      # Convert interval_num to integer and extract the corresponding row
      interval_num_int <- as.integer(interval_num)
      ranks_vector <- as.numeric(permuted_ranks[intervalNumber == interval_num_int, .SD, .SDcols = patterns("^rank_")])
      #print(ranks_vector)
      #message("Checkpoint 4 inside permuted_rank_products")
      
      # Calculate the geometric mean of ranks
      geometric_mean_rank(ranks_vector)
    })
    
    return(permuted_rank_products)
  })
  
  # Convert the results list to a tibble
  results_tibble <- tibble::tibble(interval = data$intervalNumber, 
                                   permuted_rank_products = results_list) %>% 
    unnest_wider(permuted_rank_products, names_sep = "_")
  
  print("Results tibble:")  # Debug print
  print(results_tibble)
  
  return(results_tibble)
}

# Example usage
nonZeroLocations <- locationBackbone %>%
  dplyr::select(intervalNumber, intervalChrom, intervalStart, intervalEnd, interval_rank_product, geo_mean_rank, gene_count, dplyr::starts_with("rank_")) %>%
  dplyr::filter(geo_mean_rank != length(locationBackbone$intervalNumber)) %>%
  dplyr::mutate(across(starts_with("rank_"), as.numeric))

#print("Non-zero locations:")  # Debug print
#print(nonZeroLocations)
headNonZero <- head(nonZeroLocations, n = 1000)
nPerm <- 2
permuted_rank_products_tibble <- calculate_permuted_rank_products(headNonZero, num_permutations = nPerm)

print( Sys.time() - start )

###############################
### OPTIMIZATION VERSION ### V3!!!!
###############################
start <- Sys.time()


if (.Platform$OS.type == "unix"){
  plan(multicore)
} else if  (.Platform$OS.type == "windows") {
  plan(multisession)
} else {
  plan(sequential)
}

calculate_permuted_rank_products <- function(data, num_permutations = 1000) {
  # Convert input data to data.table once for performance
  setDT(data)
  
  # Helper function to calculate geometric mean of ranks
  geometric_mean_rank <- function(ranks) {
    ranks <- as.numeric(ranks)  # Explicitly convert to numeric
    if (any(is.na(ranks))) {
      stop("Ranks contain NA values")
    }
    exp(mean(log(ranks)))  # Calculate geometric mean of ranks
  }
  
  # Iterate over each intervalNumber in the data and perform permutations sequentially
  results_list <- map(data$intervalNumber, function(interval_num) {
    # Filter current interval and same gene count intervals
    current_interval <- data[intervalNumber == interval_num, ]
    same_gene_count_intervals <- data[gene_count == current_interval$gene_count]
    
    # Perform permutation for each interval
    permuted_rank_products <- replicate(num_permutations, {
      # Sample the ranks within the groups of intervalNumber
      permuted_ranks <- same_gene_count_intervals[, lapply(.SD, function(x) x[sample(.N)]), 
                                                  by = intervalNumber, .SDcols = patterns("^rank_")]
      
      # Extract the corresponding row based on interval_num and calculate the geometric mean
      ranks_vector <- as.numeric(permuted_ranks[intervalNumber == interval_num, .SD, .SDcols = patterns("^rank_")])
      geometric_mean_rank(ranks_vector)
    })
    
    return(permuted_rank_products)
  })
  
  # Convert the results list to a tibble
  results_tibble <- tibble(interval = data$intervalNumber, 
                           permuted_rank_products = results_list) %>% 
    unnest_wider(permuted_rank_products, names_sep = "_")
  
  return(results_tibble)
}

# Example usage
nonZeroLocations <- locationBackbone %>%
  dplyr::select(intervalNumber, intervalChrom, intervalStart, intervalEnd, interval_rank_product, geo_mean_rank, gene_count, dplyr::starts_with("rank_")) %>%
  dplyr::filter(geo_mean_rank != length(locationBackbone$intervalNumber)) %>%
  dplyr::mutate(across(starts_with("rank_"), as.numeric))

#print("Non-zero locations:")  # Debug print
#print(nonZeroLocations)
headNonZero <- head(nonZeroLocations, n = 1000)
nPerm <- 2
permuted_rank_products_tibble <- calculate_permuted_rank_products(headNonZero, num_permutations = nPerm)

print( Sys.time() - start )


###############################
### OPTIMIZATION VERSION ### V4!!!!
###############################
start <- Sys.time()


if (.Platform$OS.type == "unix"){
  plan(multicore)
} else if  (.Platform$OS.type == "windows") {
  plan(multisession)
} else {
  plan(sequential)
}

calculate_permuted_rank_products <- function(data, num_permutations = 1000) {
  setDT(data)
  
  geometric_mean_rank <- function(ranks) {
    ranks <- as.numeric(ranks)
    if (any(is.na(ranks))) stop("Ranks contain NA values")
    exp(mean(log(ranks)))
  }
  
  results_list <- future_map(data$intervalNumber, function(interval_num) {
    current_interval <- data[intervalNumber == interval_num, ]
    same_gene_count_intervals <- data[gene_count == current_interval$gene_count]
    
    permuted_rank_products <- replicate(num_permutations, {
      permuted_ranks <- same_gene_count_intervals[, lapply(.SD, function(x) x[sample(.N)]), 
                                                  by = intervalNumber, .SDcols = patterns("^rank_")]
      
      ranks_vector <- as.numeric(permuted_ranks[intervalNumber == interval_num, .SD, .SDcols = patterns("^rank_")])
      geometric_mean_rank(ranks_vector)
    })
    
    return(permuted_rank_products)
  }, .options = furrr_options(seed = TRUE))  # Ensure reproducibility with seed
  
  results_tibble <- tibble(interval = data$intervalNumber, 
                           permuted_rank_products = results_list) %>% 
    unnest_wider(permuted_rank_products, names_sep = "_")
  
  return(results_tibble)
}

# Example usage
nonZeroLocations <- locationBackbone %>%
  dplyr::select(intervalNumber, intervalChrom, intervalStart, intervalEnd, interval_rank_product, geo_mean_rank, gene_count, dplyr::starts_with("rank_")) %>%
  dplyr::filter(geo_mean_rank != length(locationBackbone$intervalNumber)) %>%
  dplyr::mutate(across(starts_with("rank_"), as.numeric))

#print("Non-zero locations:")  # Debug print
#print(nonZeroLocations)
headNonZero <- head(nonZeroLocations, n = 1000)
nPerm <- 2
permuted_rank_products_tibble <- calculate_permuted_rank_products(headNonZero, num_permutations = nPerm)

print( Sys.time() - start )

#'[STATISTICS]
#Add row that counts how many values in the permuted RP tibble are <= the real RP value (RPperm <= RPreal)
#We will call this the cValue
headNonZero2 <- headNonZero %>%
  rowwise() %>%
  mutate(
    cValue = sum(
      permuted_rank_products_tibble %>%
        dplyr::filter(intervalNumber == interval) %>%
        dplyr::select(starts_with("permuted_rank_products_")) %>%
        unlist() <= interval_rank_product
    )
  ) %>%
  ungroup()

#Add row that calculates the expected RP value (Erp = cValue/number of permutations)
#Add row that calculates the p-value, in this case the percentage of false positives (PFP), where PFP = Erp/rank
headNonZero2 <- headNonZero2 %>%
  dplyr::mutate(expected_RP = (cValue) / (nPerm)) %>%
  dplyr::mutate(p_value = expected_RP/geo_mean_rank) %>%
  dplyr::mutate(fdr_corrected_p_value = p.adjust(p_value, method = "fdr"))
  
  


