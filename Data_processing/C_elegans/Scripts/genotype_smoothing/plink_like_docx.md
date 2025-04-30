## Script Overview

The script is designed to clean and preprocess genetic marker data (SNPs) from a CSV file. The dataset contains rows representing genetic markers and columns representing:

- Marker metadata (Chr, Locus, cM, Mb).
- Genotype values for 199 individuals (WN1001 to WN1199), with values like 0.0, 0.5, or 1.0 (representing homozygous reference, heterozygous, or homozygous alternate alleles, respectively).

The goal is to:

- Remove low-quality or uninformative markers.
- Impute missing genotypes.
- Reduce redundancy by pruning highly correlated markers.
- Save the cleaned dataset.

### 01 Import libraries
```python 
import pandas as pd
import numpy as np
from sklearn.impute import SimpleImputer
#from scipy.stats import chi2_contingency

```
- `pandas`: loading and data manipulation 
- `numpy`: numerical operations 
- `sklearn.impute.SimpleImputer`: A tool from scikit-learn to fill in missing values in the data

### 02 Load data 
```python 
data = pd.read_csv('genotypes.csv', low_memory=False)

```
- data is load into memory to be processed (dataframe) 

### 03 Extract genotype columns 
```python 
genotype_cols = [col for col in data.columns if col.startswith('WN')]
genotypes = data[genotype_cols]

```
- Identifies all columns in the DataFrame that start with 'WN' (e.g., WN1001, WN1002, etc.), which correspond to the genotype data for each individual.
- Creates a new DataFrame genotypes containing only these genotype columns (excluding metadata like Chr, Locus, etc) 

- How it works: 
  - The list comprehension `[col for col in data.columns if col.startswith('WN')]` loops through all column names and selects those starting with 'WN'.
  - `data[genotype_cols]` extracts only those columns into a new DataFrame genotypes.
  - Example: If the original data has columns `[Chr, Locus, cM, Mb, WN1001, WN1002, ...]`, then genotype_cols will be `['WN1001', 'WN1002', ...]`, and genotypes will include only the genotype columns.
 
### 04 Quality Control: Removes markers with high missing data 
```python 
missing_rate = genotypes.isna().mean(axis=1)
data = data[missing_rate < 0.1]

```
- Calculates the proportion of missing values (NaN) for each marker (row) across all individuals.
- Removes markers where more than 10% of the genotype values are missing. 
- Markers with many missing values are unreliable because they lack sufficient data to analyze. Removing them reduces noise and improves data quality. 

- How it works: 
  - `genotypes.isna()` creates a DataFrame of True/False values (True where a value is NaN).
  - `.mean(axis=1)` computes the average number of True values per row (i.e., the fraction of missing values for each marker).
  - `missing_rate < 0.1` creates a boolean mask where True indicates markers with less than 10% missing values.
  - `data = data[missing_rate < 0.1]` keeps only the rows (markers) where the missing rate is below 10%.
  - Example: If a marker has 199 individuals and 20 are NaN, the missing rate is 20/199 ≈ 0.1. If it’s >0.1, the marker is removed.

### 05 Quality Control: Calculate MAF and remove low MAF markers 
```python 
def calculate_maf(row):
    alleles = row.value_counts(normalize=True)
    freqs = [alleles.get(x, 0) for x in [0.0, 0.5, 1.0]]
    p = (2 * freqs[0] + freqs[1]) / (2 * sum(freqs))
    q = 1 - p
    return min(p, q)

mafs = genotypes.apply(calculate_maf, axis=1)
data = data[mafs >= 0.01]

```
- Defines a function `calculate_maf` to compute the Minor Allele Frequency (MAF) for each marker.
- Removes markers with MAF less than 0.01 (1%). 
- MAF measures how common the less frequent allele is. Markers with very low MAF (e.g., <1%) have little genetic variation, making them less informative and more prone to errors. Removing them focuses on more variable markers.

- How it works: 
  - Function `calculate_maf`:

    - `row.value_counts(normalize=True)` counts the frequency of each genotype (0.0, 0.5, 1.0) in the row, normalized to sum to 1.
    - `freqs = [alleles.get(x, 0) for x in [0.0, 0.5, 1.0]]` gets the frequencies of 0.0 (homozygous reference, AA), 0.5 (heterozygous, Aa), and 1.0 (homozygous alternate, aa). If a genotype is absent, its frequency is 0.
    - `p = (2 * freqs[0] + freqs[1]) / (2 * sum(freqs))` calculates the frequency of the reference allele (A):
        - Each 0.0 genotype contributes 2 A alleles.
        - Each 0.5 genotype contributes 1 A allele.
        - The denominator 2 * sum(freqs) accounts for the total number of alleles (2 per individual).
    - `q = 1 - p` is the frequency of the alternate allele (a).
    - `min(p, q)` returns the MAF, the frequency of the less common allele.

   - `genotypes.apply(calculate_maf, axis=1)` applies the function to each row, producing a Series of MAF values.
   - `mafs >= 0.01` creates a boolean mask for markers with MAF ≥ 1%.
   - `data = data[mafs >= 0.01]` keeps only those markers.

   - Example:

    - Suppose a marker has 100 individuals: 80 with 0.0 (AA), 15 with 0.5 (Aa), 5 with 1.0 (aa).
    - Frequencies: `freqs = [0.8, 0.15, 0.05]`.
    - Allele A count: `(2 * 80 + 15) = 175`.
    - Total alleles: `2 * 100 = 200`.
    - `p = 175/200 = 0.875`, `q = 1 - 0.875 = 0.125`.
    - MAF = `min(0.875, 0.125) = 0.125` (12.5%).
    - Since 0.125 ≥ 0.01, this marker is kept.

### 06 Quality Control: Remove monomorphic markers
```python 
std = genotypes.std(axis=1)
data = data[std > 0]

``` 
- Calculates the standard deviation of genotype values for each marker.
- Removes markers where the standard deviation is 0 (i.e., all genotypes are identical).
- Monomorphic markers (where all individuals have the same genotype, e.g., all 0.0) provide no genetic variation and are useless for most analyses. Removing them reduces redundancy.

- How it works: 
  - `genotypes.std(axis=1)` computes the standard deviation of genotype values (0.0, 0.5, 1.0) across individuals for each marker.
  - If `std = 0`, all values are the same (e.g., all 0.0), indicating a monomorphic marker.
  - `std > 0` creates a boolean mask for markers with some variation.
  - `data = data[std > 0]` keeps only those markers.
- Example: 
  - Marker with genotypes `[0.0, 0.0, 0.0, ...]` has `std = 0` and is removed.
  - Marker with genotypes `[0.0, 0.5, 1.0, ...]` has `std > 0` and is kept.
  
### 07 Imputation of missing values 
```python
imputer = SimpleImputer(strategy='most_frequent')
genotypes_imputed = pd.DataFrame(imputer.fit_transform(genotypes), columns=genotype_cols)
data[genotype_cols] = genotypes_imputed

```
- Uses `SimpleImputer` to fill in any remaining missing values (`NaN`) in the genotype data.
- Replaces missing values with the most frequent genotype (mode) for each marker.
- Updates the original `data` DataFrame with the imputed genotypes.
- Missing values can disrupt analyses. Imputing them ensures a complete dataset, improving reliability. Using the most frequent genotype is a simple approach, though it may not account for genetic patterns like LD.

- How it works: 
  - `SimpleImputer(strategy='most_frequent')` creates an imputer object that replaces NaN values with the mode of each column.
  - `imputer.fit_transform(genotypes)` computes the mode for each genotype column and fills in missing values, returning a NumPy array.
  - `pd.DataFrame(..., columns=genotype_cols)` converts the array back to a DataFrame with the same column names.
  - `data[genotype_cols] = genotypes_imputed` updates the genotype columns in the original data DataFrame with the imputed values.
- Example:
  - Suppose a marker has genotypes [0.0, NaN, 0.5, 0.0, 0.0].
  - The mode is 0.0 (most frequent).
  - The imputed genotypes become [0.0, 0.0, 0.5, 0.0, 0.0].

### 08 LD Pruning (simplified correlation based) 
```python 
def ld_pruning(genotypes, threshold=0.8, window=100):
    keep = []
    for i in range(0, len(genotypes), window):
        window_data = genotypes.iloc[i:i+window]
        corr = window_data.corr().abs()
        np.fill_diagonal(corr.values, 0)
        # Greedily select markers with low correlation
        selected = []
        for j in range(len(window_data)):
            if all(corr.iloc[j, k] < threshold for k in selected):
                selected.append(j)
        keep.extend(window_data.index[selected])
    return keep

keep_indices = ld_pruning(genotypes_imputed)
data = data.loc[keep_indices]

```
- Defines a function `ld_pruning` to remove markers that are highly correlated (in high linkage disequilibrium, LD) within a window of 100 markers.
- Keeps only markers with pairwise correlations below a threshold (0.8).
- Updates the `data` DataFrame to include only the selected markers.
- Markers in high LD (highly correlated) provide redundant information because they tend to be inherited together. Pruning them reduces dataset size and noise without losing much genetic information.

- How it works: 
  - Function `ld_pruning`: 
    - `for i in range(0, len(genotypes), window)` loops over the dataset in chunks of 100 markers (e.g., rows 0–99, 100–199, etc.).
    - `window_data = genotypes.iloc[i:i+window]` extracts the genotype data for the current window.
    - `corr = window_data.corr().abs()` computes the absolute Pearson correlation matrix between markers in the window. Correlation measures how similar the genotype patterns are.
    - `np.fill_diagonal(corr.values, 0)` sets the diagonal (correlation of a marker with itself) to 0 to avoid selecting a marker based on its self-correlation.
    - Greedy selection: 
      - `selected = []` keeps track of indices of markers to retain.
      - `for j in range(len(window_data))` loops through each marker in the window.
      - `if all(corr.iloc[j, k] < threshold for k in selected)` checks if the current marker’s correlation with all previously selected markers is below 0.8.
      - If true, the marker is added to `selected`.
      - This ensures that retained markers are not highly correlated with each other.
    - `keep.extend(window_data.index[selected])` adds the indices of selected markers to the keep list.
  - `keep_indices = ld_pruning(genotypes_imputed)` runs the function on the imputed genotypes, returning a list of row indices to keep.
  - `data = data.loc[keep_indices]` keeps only the rows (markers) with those indices.

### 09 Save cleaned data 
```python 
data.to_csv('cleaned_genotypes.csv', index=False)

```

