import argparse 
from celegans_geno_converter.converter import convert_genotype_files 

def main():
    parser = argparse.ArgumentParser(description="Convert genotype and marker files for Celegans datasets into GEMMA-compatible format.")

    parser.add_argument("geno_file", help="Path to the genotype file (e.g., file_map.txt)")
    parser.add_argument("marker_file", help="Path to the marker file (e.g., file_marker.txt)")
    parser.add_argument("output_file", help="Path to save the processed output (e.g., processed/C_elegans-Genotypes.tsv)")

    args = parser.parse_args()

    convert_genotype_files(args.geno_file, args.marker_file, args.output_file)

if __name__ == "__main__":
    main()

