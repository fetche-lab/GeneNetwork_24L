### A simple python script to extract gene IDs from a phenotype dataset and use them to query information from TAIR website for the metadata

```python
# import pandas 
import pandas as pd

#Extract the rownames and save them in a separate csv file 
model_df = pd.read_csv('phenotype_im_Harm_arabid.csv')

#Inspect the dataset 
model_df.shape
model_df.head(20)
model_df.tail(20)

#extract the rownames and save them in a separate text file
model_ids = [data for data in model_df['AraSeqID']]
print(len(model_ids)) # the value should equal the one in model_df.shape, on the x axis

#specify the filename 
file_name = 'Harm_model_id.txt'

#write to the txt file 
with open(file_name, 'w') as file1:
    
    #write each row from the list 
    for data in model_ids: 
        file1.write(data + '\n')
print(f'The data has been written to {file_name}')

##Dataframe containing the downloaded information containing gene description
meta1_df = pd.read_table('Gene_description-harm.txt') 

## create a new dataframe with only the selected columns 
meta2_df = pd.DataFrame(meta1_df[['Locus Identifier', 'Gene Model Name', 'Gene Model Description','Primary Gene Symbol']])

#inspect the dataframe 
meta2_df.head()                  

```
