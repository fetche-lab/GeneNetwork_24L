#!/usr/bin/env python3 

## import packages 
import argparse
import pandas as pd

## the program 
def marker_cols_rearrange(marker_file, sheet_name, output_file):
    """ 
    Checks and rearrange the columns in marker file into a good order 
    
    Parameters: 
    - Marker_file(str): Path to the marker file (e.g marker_file.xlsx)
    - sheet_name(str): The name of your sheet in excel (e.g sheet01) 
    - Output_file(str): Path to save the processed output file  
    """
    ## load in the file 
    marker_df = pd.read_excel(marker_file, sheet_name=sheet_name)

    ## drop the first column 
    del marker_df["name"]

    ## rearrange the columns 
    marker_df1 = marker_df.set_axis(['name', 'Chromosome', 'Pos_WS258'], axis=1)

    ## save the file in txt format 
    marker_df1.to_csv(output_file, header = True, index = None, sep="\t") 

    # Exit message 
    print(f"Processing complete: File saved at: {output_file}")

def main():
    parser = argparse.ArgumentParser(description="Rearranges columns in marker file in the right order")
    parser.add_argument("marker_file", help="Path to the marker file {it should be in Excel format}")
    parser.add_argument("sheet_name", help="Name of the sheet with your data in excel")
    parser.add_argument("output_file", help="Path to save your output file")

    args = parser.parse_args()
    marker_cols_rearrange(args.marker_file, args.sheet_name, args.output_file)

if __name__ == "__main__":
    main()
