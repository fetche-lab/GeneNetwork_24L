### Processing Kilifish data for GeneNetwork 
The following contains scripts used to process Kilifish datasets 

#### 01 Genotype dataset 
```python 
## import packages 
import pandas as pd 

## load in the file into dataframe 
geno_df = pd.read_csv("killi.v0.1.cross", low_memory=False, header=None)

## examine the file to see what it contains and the general formatting 
geno_df.head()
geno_df.tail() 
geno_df.shape 
geno_df.columns 
geno_df.info() 

## in our case, this file had extra rows that needed to be removed 
### getting rid of those rows using their indeces 
geno_df1 = geno_df.drop([0, 2, 3, 4, 5])
### reset the index order for the rows 
geno_df1.reset_index(drop=True, inplace=True)
### now set the first row values to be column headers 
geno_df1.columns = geno_df1.iloc[0] #assign the 1st row as column header
geno_df1 = geno_df1[1:] #remove the first row from the data 
geno_df1.reset_index(drop=True, inplace=True) #reset the indeces once more 

## rename the column headers for Xsome and position 
geno_df2 = geno_df1.rename(columns={geno_df1.columns[1]:"Chr",  geno_df1.columns[2]:"Position"})

## inspect your df and save it to a csv file 
geno_df2.to_csv("Kilifish_genotypes.csv", index=None, header=True)

```

#### 02 Case Attributes 
```python 
## import packages 
import pandas as pd 

## load in the file into dataframe 
attr_df = pd.read_csv("2025-03-25_QTL_Master_table_visual_and_transcriptomics_proteomics_full.csv", low_memory=False)

## inspect the columns, and see what's need to be used as case attributes depending on your criteria as a scientist 
### in this case, the following columns were selected to be used as case attributes 
cols_to_copy = ['ID', 'Generation', 'DOB', 'DOS', 'sex','Liver.phenotype...1.0.', 'Split.tail...1.0.', 'Iris.colour..w.y.o.','Iris.black...1.0.', 'Upper.head.1.0', 'lower.head.1.0','Lip.protrusion.1.0', 'Gill.defect..KD..1.0.', 'body.redness..0.3.','body.turquoise..0.3.']
### create a new df with the selected columns 
attr_df1 = attr_df[cols_to_copy]
### fill missing values with the `NA` string 
attr_df1.fillna("NA", inplace=True)

## inspect the final df, then save it to a csv file 
attr_df1.to_csv("Kilifish_caseAttributes.csv", index=None, header=True)

```

#### 03 Phenotype dataset (non-expression/experimental phenotypes)

```python 
## import packages 
import pandas as pd 

## load in the file into a dataframe 
phenotype_df = pd.read_csv("2025-03-25_QTL_Master_table_visual_and_transcriptomics_proteomics_full.csv", low_memory=False)

## inspect the df, then select needed columns 
cols_to_copy = ['ID', 'Generation', 'Brain.weight', 'sex','Liver.phenotype...1.0.', 'Split.tail...1.0.','Iris.black...1.0.', 'Upper.head.1.0', 'lower.head.1.0','Lip.protrusion.1.0', 'Gill.defect..KD..1.0.', 'body.redness..0.3.','body.turquoise..0.3.', 'tail.red..', 'tail.yellow..', 'tail.black..','tail.transp..', 'CM.length.in.pixel', 'CM.2.surface.in.pixel','Fish.length.total', 'Fish.length.no.tail', 'Fish.surface.total','tail.surface.area', 'eye.surface', 'Cataract.surface', 'tail.length']

### create a new dataframe with selected columns 
phenotype_df1 = phenotype_df[cols_to_copy]
### fill in empty cells 
phenotype_df1.fillna("NA", axis=1, inplace=True)

## inspect the final df, and save it into a csv file 
phenotype_df1.to_csv("Kilifish_exp_phenotypes.csv", index=None, header=True)

```

#### 04 Phenotype dataset (log-expression data/proteomics data?)

```python 
## import packages 
import pandas as pd 

## load in data into dataframe 
log2_df = pd.read_csv("2025-03-25_QTL_Master_table_visual_and_transcriptomics_proteomics_full.csv", low_memory=False)

## inspect the df, select needed columns 
cols_to_copy = list(log2_df.columns[29:len(log2_df)]) #create a list of columns needed 
log2_df1 = log2_df[cols_to_copy] #create a new df with selected columns 
log2_df1.insert(0, "IDs", log2_df["ID"]) #insert the rownames 
log2_df1.fillna("NA", inplace=True) #fill empty cells 

## inspect the final df, save it into a csv file 
log2_df1.to_csv("Kilifish_expressionData.csv", header=True, index=None)

```
#### next is the annotation file.., 