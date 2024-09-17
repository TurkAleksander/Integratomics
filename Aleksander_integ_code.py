import os
import pandas as pd

###################################
########## SETUP ##################
###################################
# Set workdir, intialize study info dataframe
# Read the headers of the files to extract study information

directory = "C:/Users/gkgen69/Desktop/VS_code/Integ_MS/Code/TestData"
setwkdir = os.chdir(directory)
studyInfoDF = pd.DataFrame(columns=["fileName", "Study_type", "Study_name", "Data_type"])

# Extract study information from the first three lines of each file
for filename in os.listdir(directory):
    if filename.endswith(".txt"):
        file_path = os.path.join(directory, filename)
        
        # Read the first three lines of the file
        with open(file_path, 'r') as file:
            lines = file.readlines()
            line1 = lines[0].strip() if len(lines) > 0 else ""
            line2 = lines[1].strip() if len(lines) > 0 else ""
            line3 = lines[2].strip() if len(lines) > 0 else ""
        
        # Append the data to the DataFrame
        studyInfoDF = pd.concat([studyInfoDF, pd.DataFrame([{
            "fileName": filename,
            "Study_type": line1,
            "Study_name": line2,
            "Data_type": line3
            }])], ignore_index=True)	

# Display the DataFrame
print(studyInfoDF)

#Import base file with chromosome lengths, data acquired from https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&chromInfoPage=
hg38BaseFile = pd.read_csv("C:/Users/gkgen69/Desktop/VS_code/Integ_MS/Code/hg38_UCSC_chrom_lengths.txt", sep="\t", header = None)

###################################
########## MAKE BACKBONE ##########
###################################
#Note: the backbone in this case refers to the genomic locations that will be used for summing the signals across various intervals
#For its construction it requires data on chromosome lengths, which is stored in the hg38BaseFile


#If it already exists, read it, if not, make it
if (os.path.exists("C:/Users/gkgen69/Desktop/VS_code/Integ_MS/Code/locationBackbone.txt") == False):
    locationBackbone = pd.DataFrame(columns=["intervalNumber","intervalChrom", "intervalStart", "intervalEnd"])
    #make backbone for genomic locations
    #hg38BaseFile = hg38BaseFile[:1]

    intervals = []

    # Initialize intervalCounter and other parameters
    intervalCounter = 0
    step = 10000
    overlap = 5000

    # Iterate through hg38BaseFile to generate intervals
    for entry in hg38BaseFile.iterrows():
        chrom = entry[1][0]
        length = entry[1][1]
        print(entry)
        for i in range(1, length, step - overlap):
            intervalCounter += 1

            intervalStart = i
            intervalEnd = i + step - 1

            # Collect the interval in the list
            intervals.append({
                "intervalNumber": intervalCounter,
                "intervalChrom": chrom,
                "intervalStart": intervalStart,
                "intervalEnd": intervalEnd
            })

    # Create the locationBackbone DataFrame from the collected intervals
    locationBackbone = pd.DataFrame(intervals, columns=["intervalNumber", "intervalChrom", "intervalStart", "intervalEnd"])

    # Print the first 10 rows of the DataFrame and save the DF to a file
    print(locationBackbone[:10])
    locationBackbone.to_csv("C:/Users/gkgen69/Desktop/VS_code/Integ_MS/Code/locationBackbone.txt", sep="\t", index=False)
else:
    locationBackbone = pd.read_csv("C:/Users/gkgen69/Desktop/VS_code/Integ_MS/Code/locationBackbone.txt", sep="\t", header = 0)

###################################
### ADD STUDIES TO BACKBONE #######
###################################

#locationBackboneStudies = locationBackbone.copy()

#Try it on a smaller subset of cata
locationBackboneStudies = locationBackbone[locationBackbone["intervalChrom"] == "chr6"]
print(locationBackboneStudies[:10])

# Iterate through the unique study types in studyInfoDF
for studyType in studyInfoDF["Study_type"].unique():
    print(studyType)
    # Get the files for the current study type
    studyFiles = studyInfoDF[studyInfoDF["Study_type"] == studyType]["fileName"]
    print(studyFiles)
    
    # Initialize an empty DataFrame to hold all study data for the current study type
    all_study_data = pd.DataFrame(columns=['chrom', 'start', 'end', 'value'])
    
    # Iterate through the files to read them and drop the first three lines (header)
    for file in studyFiles:
        print(file)
        with open(file, 'r') as f:
            lines = f.readlines()
            lines = lines[3:]
        
        # Assuming the study data is tab-separated and has columns: chrom, start, end, value
        study_data = pd.DataFrame([line.strip().split('\t') for line in lines], columns=['chrom', 'start', 'end', 'value'])
        study_data['start'] = pd.to_numeric(study_data['start'], errors='coerce').fillna(0).astype(int)
        study_data['end'] = pd.to_numeric(study_data['end'], errors='coerce').fillna(0).astype(int)
        study_data['value'] = pd.to_numeric(study_data['value'], errors='coerce').fillna(0.0).astype(float)
        
        # Append the current study data to the all_study_data DataFrame
        all_study_data = pd.concat([all_study_data, study_data], ignore_index=True)
    
    # Initialize the new column with 0.0 (float)
    locationBackboneStudies[f'{studyType}_value'] = 0.0
    
    # Use vectorized operations to sum the values for each interval
    for index, row in locationBackboneStudies.iterrows():
        chrom = row['intervalChrom']
        start = row['intervalStart']
        end = row['intervalEnd']
        
        # Find matching study data
        matching_data = all_study_data[(all_study_data['chrom'] == chrom) & 
                                       ((all_study_data['start'] <= end) & 
                                       (all_study_data['end'] >= start)
                                       |
                                       (all_study_data['start'] >= start) &
                                       (all_study_data['start'] <= end)
                                       |
                                       (all_study_data['end'] >= start) &
                                       (all_study_data['end'] <= end))]            
        # If matching data is found, sum the values and attach to a new column
        if not matching_data.empty:
            summed_value = matching_data['value'].sum()
            locationBackboneStudies.loc[index, f'{studyType}_value'] += summed_value


print(locationBackboneStudies.head())
print(locationBackboneStudies[:10])
locationBackboneStudies.to_csv("C:/Users/gkgen69/Desktop/VS_code/Integ_MS/Code/locationBackboneTest.txt", sep="\t", index=False)

#Old version:
for studyType in studyInfoDF["Study_type"].unique():
    print(studyType)
    # Get the files for the current study type
    studyFiles = studyInfoDF[studyInfoDF["Study_type"] == studyType]["fileName"]
    print(studyFiles)
    
    # Iterate through the files to read them and drop the first three lines (header)
    for file in studyFiles:
        print(file)
        with open(file, 'r') as f:
            lines = f.readlines()
            lines = lines[3:]
        
        # Assuming the study data is tab-separated and has columns: chrom, start, end, value
        study_data = pd.DataFrame([line.strip().split('\t') for line in lines], columns=['chrom', 'start', 'end', 'value'])
        study_data['start'] = pd.to_numeric(study_data['start'], errors='coerce').fillna(0).astype(int)
        study_data['end'] = pd.to_numeric(study_data['end'], errors='coerce').fillna(0).astype(int)
        study_data['value'] = pd.to_numeric(study_data['value'], errors='coerce').fillna(0.0).astype(float)
        
        # Initialize the new column with 0
        locationBackboneStudies[f'{studyType}_value'] = 0.0
        
        # Merge the study data with locationBackboneStudies based on genomic locations
        for index, row in locationBackboneStudies.iterrows():
            chrom = row['intervalChrom']
            start = row['intervalStart']
            end = row['intervalEnd']
            
            # Find matching study data
            #matching_data = study_data[(study_data['chrom'] == chrom) & 
            #                           (study_data['start'] <= end) & 
            #                           (study_data['end'] >= start)]
            
            # Find matching study data - the chromosome has to match AND
            #                          - the study data has to be within the backbone interval
            #                          - OR
            #                          - the start or end of the study data has to be within the backbone interval
            # This is done to account for the fact that backbone intervals are 5000 bp in length and the study data can be longer
            matching_data = study_data[(study_data['chrom'] == chrom) & 
                                       ((study_data['start'] <= end) & 
                                       (study_data['end'] >= start)
                                       |
                                       (study_data['start'] >= start) &
                                       (study_data['start'] <= end)
                                       |
                                       (study_data['end'] >= start) &
                                       (study_data['end'] <= end))]            
            # If matching data is found, attach the value to a new column
            if not matching_data.empty:
                summed_value = matching_data['value'].sum()
                locationBackboneStudies.loc[index, f'{studyType}_value'] += summed_value

print(matching_data.head())
print(locationBackboneStudies[:10])
locationBackboneStudies.to_csv("C:/Users/gkgen69/Desktop/VS_code/Integ_MS/Code/locationBackboneTest.txt", sep="\t", index=False)

print(studyType)

#I think Aleš's method isn't the same as Breitling's? Check it in detail.

#Iterate through a list of studies with the same study type

#Detect data type

#Extract data from the files

#Detect p-values or already -log10(p), transform if necessary

#Sum across each interval

#Do this for each study type

#Make a new column in a separate dataframe with interval ranks per study type

#Calculate rank product per interval across all groups (aleš, eq. 1)