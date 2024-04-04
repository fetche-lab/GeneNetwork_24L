'''
The following python script includes python codes that were used for data wrangling of medaka fish genotype file. 

Processes carried out included: 
01 Adding sample ids as column headers 
02 Adding marker ids as rownames/ indeces 
03 Converting the genotype values into the gemma encoded format 
04 Divide the dataset into chunks for memory efficiency while doing the analysis

'''
## read the information in the file with headers, then extract them as a list
headers_lst = []
with open('mikk_headers.txt', 'r') as file1: 
    headers = file1.readlines()
for line in headers: 
    headers_lst.append(line[:-1]) # this removes the last `\n` character 
print(len(headers_lst)) 

## Adding the genotype markers created from pos+chrom 
marker_ids = []
with open('./mikk_clean_markers_ed2.txt', 'r') as file:
    markers = file.readlines()
for lines in markers: 
    marker_ids.append(lines[:-1]) 
print(marker_ids[0:10])

### Updated version of the code 
import pandas as pd

# Define the file path and column names (if applicable)
file_path = './mikk_clean_filtrd_genotypes.txt' 
chunk_size = 10**5  # Adjust this according to your system's memory capacity

# Initialize an empty list to store processed chunks
processed_chunks = []

# Read the file in chunks
for chunk in pd.read_csv(file_path, header = None, sep = '\t', chunksize=chunk_size):
    
    # Step 1: Remove the last column
    chunk.drop(columns=[chunk.columns[-1]], inplace=True)
    
    # Step 2: Add new column headers
    chunk.columns = headers_lst # Replace with your new column headers
    
    # Step 3: Replace some values in the dataset
    chunk.replace({'0/3':'0.0','2/6':'1.0', '1/2':'1', '1/4':'1', '3/3':'1.0',  
                  '4/4':'1.0', '0/4':'0.0', '0/5':'0.0', '2/4':'1.0', '6/6':'1.0', 
                  '5/5':'1.0', '2/3':'0.0', '0/6':'0.0', '4/6':'1.0', '3/4':'1.0', 
                  '2/5':'1.0', '1/6':'1.0', '3/5':'1.0', '1/5':'1.0', '3/6':'1.0', 
                  '1/1':'1', '0/0':'0', '2/2':'1', '1/0':'0.5', '0/2':'0.0', '4/5':'1.0', 
                  '2/0':'1', '1/3':'1.0', '5/6':'1.0', '0/1':'0.5', './.':'NA'}, inplace=True) # Replace with your column name and values
    
    # Step 4: Add new row names to the datasets (optional, depending on what you mean by "rownames")
    # You can set an index if you want specific row names
    # chunk.index = ['row_name_1', 'row_name_2', ...] 
    
    # Append the processed chunk to the list
    processed_chunks.append(chunk)

# Merge the processed chunks into a final dataframe
final_df = pd.concat(processed_chunks)

# Optionally, reset the index of the final dataframe
final_df.reset_index(drop=True, inplace=True)

## adding marker ids to the final processed genotype dataframe, and inspecting the dataframe 
final_df.index = marker_ids
final_df.head()

## saving the file 
final_df.to_csv('mikk_clean_geno_marked2.csv', header=True)