## import pandas 
import pandas as pd 

## import the geno file as a dataframe 
mikk_df = pd.read_csv('./mikk_small_filtrd_genotypes.txt', sep='\t') 

## read the information in the file with headers, then extract them as a list
headers_lst = []
with open('mikk_headers.txt', 'r') as file1: 
    headers = file1.readlines()
    
for line in headers: 
    headers_lst.append(line[:-1]) # this removes the last `\n` character 
print(len(headers_lst))

## remove the column with empty values 
mikk_df1 = mikk_df.drop(columns=['Unnamed: 79'])

## append the headers corresponding to the genotype file 
mikk_df1.columns = headers_lst

## check for the unique values 
unique_lst=[]
for names in mikk_df1.columns: 
    unique_names = mikk_df1[names].unique()
    unique_lst.append(unique_names)
for data in unique_lst:
    print(data)
    
## replace values in the geno file with the gemma required standard
#df.replace({'old_value1': 'new_value1', 'old_value2': 'new_value2'}, inplace=True)
mikk_df1.replace({'0/3':'0.0','2/6': '1.0', '1/2':'1', '1/4':'1', '3/3':'1.0',  
                  '4/4':'1.0', '0/4':'0.0', '0/5':'0.0', '2/4':'1.0', '6/6':'1.0', 
                  '5/5':'1.0', '2/3':'0.0', '0/6':'0.0', '4/6':'1.0', '3/4':'1.0', 
                  '2/5':'1.0', '1/6':'1.0', '3/5':'1.0', '1/5':'1.0', '3/6':'1.0', 
                  '4/5':'1.0', '2/0':'1', '1/3':'1.0', '0/0':'0', '0/1': '0.5',
                  '1/1':'1', '0/1':'0.5', '1/2':'1','5/6':'1.0'}, inplace=True)

## save the dataframe into a csv file format 
mikk_df1.to_csv('mikk_small_genotype_gemma_encoded.csv', index = None, header = True)


