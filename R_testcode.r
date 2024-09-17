install.packages("dplyr")
install.packages("ggplot2")
install.packages("readr")
install.packages("stringr")
install.packages("tibble")

library(dplyr)
library(ggplot2)
library(readr)
library(stringr)
library(tibble)

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
  locationBackbone <- read.table("/cmg1scratch/PROJECTS/MS_integratomics/Data_integration/locationBackbone.txt", sep="\t")
}


#'[Extract study data by study type]
#

#Temporarily make a subset
studyInfoDF <- studyInfoDF %>%
  dplyr::filter(studyType == "TxBlood")

for (StudyType in unique(studyInfoDF$studyType)) {
  fileNames <- studyInfoDF %>%
    dplyr::filter(studyType == StudyType) %>%
    dplyr::select(studyFile) %>%
    unlist()
  
  print(StudyType)
  print(fileNames)
  
  # Store all signals from the same study type into one dataframe
  fileDataHolder <- tibble::tibble()
  for (file in fileNames) {
    currentFileData <- readr::read_tsv(file, skip = 3, col_names = FALSE)
    print(file)
    fileDataHolder <- bind_rows(fileDataHolder, currentFileData)
  }
  colnames(fileDataHolder) <- c("intervalChrom", "intervalStart", "intervalEnd", "intervalSignal")
  
  # Add a new column for the current StudyType to store the summed interval signals
  studyTypeColumn <- paste0("intervalSignalSum_", StudyType)
  locationBackbone <- locationBackbone %>% mutate(!!studyTypeColumn := 0)
  
  # Check and sum interval signals
  for (i in 1:nrow(fileDataHolder)) {
    chrom <- fileDataHolder$intervalChrom[i]
    start <- fileDataHolder$intervalStart[i]
    end <- fileDataHolder$intervalEnd[i]
    signal <- fileDataHolder$intervalSignal[i]
    
    locationBackbone <- locationBackbone %>%
      mutate(!!studyTypeColumn := if_else(
        intervalChrom == chrom & 
        ((intervalStart <= start & intervalEnd >= start) | 
         (intervalStart <= end & intervalEnd >= end)),
        !!sym(studyTypeColumn) + signal,
        !!sym(studyTypeColumn)
      ))
  }
}

#Copilot suggestion:
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
  
  # Add a new column for the current StudyType to store the summed interval signals
  studyTypeColumn <- paste0("intervalSignalSum_", StudyType)
  locationBackbone[, (studyTypeColumn) := 0]
  
  # Perform the join and sum the signals
  fileDataHolder[, intervalSignal := as.numeric(intervalSignal)]
  summedSignals <- fileDataHolder[locationBackbone, on = .(intervalChrom = intervalChrom), 
                                  nomatch = 0, allow.cartesian = TRUE][
    intervalStart <= i.intervalEnd & intervalEnd >= i.intervalStart, 
    .(intervalSignalSum = sum(intervalSignal)), 
    by = .EACHI]
  
  # Update the locationBackbone with the summed signals
  locationBackbone[summedSignals, on = .(intervalChrom, intervalStart, intervalEnd), 
                   (studyTypeColumn) := intervalSignalSum]
}

#'[ADD GEOMETRIC MEAN, RANKS]
# Identify columns whose name contains "signal", apply geometric mean across intervals
signal_cols <- grep("Signal", names(locationBackbone), value = TRUE)
geo_mean <- function(x) {
  exp(mean(log(x)))
}
locationBackbone[, geo_mean_signal := apply(.SD, 1, geo_mean), .SDcols = signal_cols]


# Add a ranking column where ties are assigned the highest rank
# Example : if intervals A, B and C are tied for 1st place, they will each be assigned the rank "3"
locationBackbone[, geo_mean_rank := frank(-geo_mean_signal, ties.method = "max")]


#'[ADD PERMUTATIONS]

#Find the locations of intervals with the same gene density, for the data on gene density we will use the Ensembl GRCh38.p14 assembly
#We included protein-coding genes from chr 1-22+X+Y+MT, that have a UCSC Gene Stable ID (20 037 genes)
biomartGeneLocations <- read.table("/cmg1scratch/PROJECTS/MS_integratomics/Data_integration/mart_export.txt", sep="\t", header = TRUE) %>%
  dplyr::distinct() %>%
  dplyr::select(Chromosome.scaffold.name, Gene.start..bp., Gene.end..bp., Gene.name, Gene.stable.ID) %>%
  dplyr::rename(Chrom = Chromosome.scaffold.name, Start = Gene.start..bp., End = Gene.end..bp., Gene_name = Gene.name, Ensembl_ID = Gene.stable.ID) %>%
  data.table::as.data.table()

biomartGeneLocations[, Chrom := paste0("chr", Chrom)]


setDT(biomartGeneLocations)

# Set keys for interval join
setkey(locationBackbone, intervalChrom, intervalStart, intervalEnd)
setkey(biomartGeneLocations, Chrom, Start, End)

# Perform the interval join and count the number of genes in each interval
overlaps <- foverlaps(biomartGeneLocations, locationBackbone, by.x = c("Chrom", "Start", "End"), by.y = c("intervalChrom", "intervalStart", "intervalEnd"), type = "any", nomatch = 0)
geneCounts <- overlaps[, .N, by = .(intervalChrom, intervalStart, intervalEnd)]

# Add the count of genes to locationBackbone
locationBackbone[, gene_count := 0]
locationBackbone[geneCounts, on = .(intervalChrom, intervalStart, intervalEnd), gene_count := N]


geo_mean_rank <- locationBackbone$geo_mean_rank
unique_values <- unique(geo_mean_rank)
n_positions <- length(geo_mean_rank)
n_unique_values <- length(unique_values)

# Function to check if 99% of positions have seen all unique values
positions_covered <- function(seen_matrix) {
  covered_positions <- apply(seen_matrix, 1, function(x) all(x))
  mean(covered_positions) >= 0.99
}

# Monte Carlo simulation to estimate the number of shuffles needed
estimate_shuffles <- function(tolerance = 0.05) {
  shuffle_counts <- numeric()
  sim <- 0
  
  while (TRUE) {
    sim <- sim + 1
    seen_matrix <- matrix(FALSE, nrow = n_positions, ncol = n_unique_values)
    shuffle_count <- 0
    print(paste("Simulation:", sim))
    
    while (!positions_covered(seen_matrix)) {
      shuffle_count <- shuffle_count + 1
      shuffled_geo_mean_rank <- sample(geo_mean_rank)
      
      for (pos in 1:n_positions) {
        value <- shuffled_geo_mean_rank[pos]
        value_index <- which(unique_values == value)
        seen_matrix[pos, value_index] <- TRUE
      }
      
      # Debugging: Print the seen_matrix occasionally
      if (shuffle_count %% 100 == 0) {
        print(paste("Shuffle count:", shuffle_count))
        print(seen_matrix)
      }
    }
    
    shuffle_counts <- c(shuffle_counts, shuffle_count)
    
    # Calculate the mean and standard error
    mean_shuffles <- mean(shuffle_counts)
    std_error <- sd(shuffle_counts) / sqrt(sim)
    
    # Debugging: Print the current mean and standard error
    print(paste("Mean shuffles:", mean_shuffles))
    print(paste("Standard error:", std_error))
    
    # Check if the standard error is within the tolerance
    if (!is.na(std_error) && std_error < tolerance) {
      break
    }
  }
  
  list(mean_shuffles = mean_shuffles, simulations = sim)
}

# Run the Monte Carlo simulation
result <- estimate_shuffles(tolerance = 0.05)
print(paste("Estimated number of shuffles needed:", result$mean_shuffles))
print(paste("Number of simulations run:", result$simulations))





geo_mean_rank <- locationBackbone$geo_mean_rank
unique_values <- unique(geo_mean_rank)
n_positions <- length(geo_mean_rank)
n_unique_values <- length(unique_values)

# Debugging: Print the initial values
print(paste("n_positions:", n_positions))
print(paste("n_unique_values:", n_unique_values))

# Monte Carlo parameters
target_percentage <- 0.99  # We want 99% coverage of all positions
tolerance <- 0.05  # Tolerance for the standard error

# Function to simulate a single Monte Carlo run
monte_carlo_run <- function() {
  seen_values <- vector("list", n_positions)  # Track values seen at each position
  for (i in 1:n_positions) {
    seen_values[[i]] <- c()  # Initialize empty list for each position
  }
  
  permutations_needed <- 0
  
  while (TRUE) {
    # Shuffle the geo_mean_rank
    shuffled_geo_mean_rank <- sample(geo_mean_rank)
    permutations_needed <- permutations_needed + 1
    
    # Debugging: Print the length of shuffled_geo_mean_rank
    print(paste("Length of shuffled_geo_mean_rank:", length(shuffled_geo_mean_rank)))
    
    # Add the shuffled values to the corresponding position trackers
    for (i in 1:n_positions) {
      # Debugging: Print the current index and value
      if (i > length(shuffled_geo_mean_rank)) {
        print(paste("Index i out of bounds:", i))
        next
      }
      print(paste("Index i:", i))
      print(paste("Value at shuffled_geo_mean_rank[i]:", shuffled_geo_mean_rank[i]))
      
      seen_values[[i]] <- unique(c(seen_values[[i]], shuffled_geo_mean_rank[i]))
    }
    
    # Check how many positions have seen all unique values
    covered_positions <- sum(sapply(seen_values, function(x) length(x) == n_unique_values))
    
    # Check if at least 99% of the positions have seen all unique values
    if (covered_positions / n_positions >= target_percentage) {
      return(permutations_needed)
    }
  }
}

# Monte Carlo simulation to estimate the number of shuffles needed
estimate_shuffles <- function(tolerance = 0.05) {
  shuffle_counts <- numeric()
  sim <- 0
  
  while (TRUE) {
    sim <- sim + 1
    shuffle_count <- monte_carlo_run()
    shuffle_counts <- c(shuffle_counts, shuffle_count)
    
    # Calculate the mean and standard error
    mean_shuffles <- mean(shuffle_counts)
    std_error <- sd(shuffle_counts) / sqrt(sim)
    
    # Print the current mean and standard error
    print(paste("Simulation:", sim))
    print(paste("Mean shuffles:", mean_shuffles))
    print(paste("Standard error:", std_error))
    
    # Check if the standard error is within the tolerance
    if (!is.na(std_error) && std_error < tolerance) {
      break
    }
  }
  
  return(mean_shuffles)
}

# Run the Monte Carlo simulation
estimated_shuffles <- estimate_shuffles(tolerance)
print(paste("Estimated number of shuffles needed:", estimated_shuffles))





