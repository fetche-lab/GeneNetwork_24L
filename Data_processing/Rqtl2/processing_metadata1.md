### Retrieving metadata information for Arabidopsis thaliana Hartanto et al dataset from gff (general feature format) file  
#### The following contains a python script which is used to parse through a gff file, and retrieve relatable descriptive information for the Hartanto et al arabid dataset from tair.org website
```python
#load the required package(s)
import pandas as pd

#load the arabid dataset with the required gene ids
data_df1 = pd.read_csv('../data/Harms_sup_data/rqlt2_data-H/phenotype_im_Harm_arabid.csv')

#inspect the dataset
print(data_df1.shape) #(29913, 46)
print(data_df1.head(20)) 
print(data_df1.iloc[:,0]) #get rownames for the first column 

#functions for parsing process 
def parse_gff_line(line):
    '''
    Extracts information about genes from the gff files by going through
    each line in the file 
    '''
    #split the line into fields using TAB as the delimiter 
    fields = line.strip().split('\t')
    
    #assign individual fields to variables 
    chrom, source, feature_type, start, end, _, strand, _, attributes = fields 
    
    #split attributes and handle cases where '=' appears in values 
    
    #create a dictionary to store attribute key-value pairs 
    attributes_dict = {} 
    
    #split the 'attributes' field using ';' as the delimiter
    for attr in attributes.split(';'):
        
        #split only on the first occurance 
        key_value = attr.split('=', 1) 
        
        #check if the splits results in two parts before adding the values to the dict 
        if len(key_value) == 2: 
            attributes_dict[key_value[0]] = key_value[1]
        
    gene_id = attributes_dict.get('ID')
    gene_name = attributes_dict.get('symbol')
    description = attributes_dict.get('Note', 'NA') 
    
    #return gene_id, gene_name, f'{chrom}:{start}-{end}({strand})', description 
    return gene_id, gene_name, chrom, f'{start}-{end}', strand, description

def parse_gff(filepath):
    ''' reads and parses the provided gff file '''
    gene_data =[]
    with open(filepath, 'rt', encoding='latin-1') as file: # rt is the text reading mode 
        for line in file:
            if line.startswith('#') or not line.strip():
                continue #skips header or empty lines 
            if '\tgene\t' in line: #process only lines with gene information 
                #yield parse_gff_line(line)
                gene_data.append(parse_gff_line(line))
                
    # Convert the list of parsed gene information into a DataFrame
    columns = ['Gene_ID', 'Gene_Name', 'Chrom', 'Start_End', 'Gene_strand','Description']
    df = pd.DataFrame(gene_data, columns=columns)
    return df

#create a dataframe from the gff file using the parse_gff function 
filepath = '../data/Harms_sup_data/Araport11_GFF3_genes_transposons.current.gff/Araport11_GFF3_genes_transposons.Jan2024.gff'
results_df = parse_gff(filepath)

#inspect the newly created dataframe
print(results_df.shape) #(33243, 6)
print(results_df.head(20))

#create another dataframe with the descriptive information corresponding to the arabidopsis dataset under study
new_metadata_df = results_df[results_df.iloc[:,0].isin(data_df1.iloc[:,0])]

#inspect the dataset
print(new_metadata_df.shape) #(26546, 6)
print(new_metadata_df.head(20))


```
