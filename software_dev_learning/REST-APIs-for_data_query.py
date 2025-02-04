'''
The following contains simple python scripts to help with consuming the rest api endpoints in ncbi, using the eutilities 
software available. 

The key features explored so far include; 
 - information about the available databases 
 - information found within a specific database 
 - the esearch tool
 - the einfo tool 
 - the esummary tool 
 - the efetch tool 

more experimentation continues as the scripter gets to decrypt the power of rest apis. This will involve integrating these 
options and the remaining ones so as to manage retrieving fully annotated dataset of interest. 
'''

## database names 
import requests 

base_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/einfo.fcgi'

params = {'retmode':'json'}

response = requests.get(base_url, params=params) 
print(response.json())



## specify the database 
import requests 
base_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/einfo.fcgi'

## parameter(s)
params = {'db':'gds',
          'retmode': 'json'
         } 

response = requests.get(base_url, params=params) 

if response.status_code == 200:
    print(response.json())
else:
    print(f'Error: {response.status_code}')
    

    
## esearch 
import requests 
base_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi'

# Specify the parameters for the API request
organism = 'Homo sapiens' 
max_results = 10 
dataset_id = 200255729
params = {
        'db': 'gds',
        'term': f'{organism}[Organism]',
        'term':'cancer',
        'retmax': max_results,
        'usehistory': 'y',
        'retmode': 'json'
    }

response = requests.get(base_url, params=params)

if response.status_code == 200:
    print(response.json())
else:
    print(f'Error: {response.status_code}')

    
## esummary
import requests 
base_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi'

# Specify the parameters for the API request
organism = 'Homo sapiens' 
max_results = 10 
dataset_id =  ['200255729', '200255707', '200255706', '200255672', '200255659', '200249500', 
               '200244327', '200244326', '200244287', '200242717']
params = {
        'db': 'gds',
        'webenv': 'MCID_65d39d215a485a235073b5a8',
        'querykey': '1',
        'id': dataset_id,
        'usehistory': 'y',
        'retmode':'json'
    }

response = requests.get(base_url, params=params)

if response.status_code == 200:
    print(response.json())
else:
    print(f'Error: {response.status_code}')

    
    
## Efetch 
import requests 
base_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'

params = {'db':'gds',
          'webenv': 'MCID_65d39d215a485a235073b5a8',
          'querykey': '1',
          'id': ['200255729', '200255707', '200255706', '200255672', 
                 '200255659', '200249500', '200244327', '200244326', 
                 '200244287', '200242717'], 
          'usehistory': 'y'}

response = requests.get(base_url, params=params)

if response.status_code == 200: 
    print(response.text)
else:
    print(f'Error: {response.status_code}')