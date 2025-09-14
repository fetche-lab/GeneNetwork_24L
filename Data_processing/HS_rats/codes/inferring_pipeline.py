#!/usr/bin/env python3 

"""
Pipeline script to: 
    1. Split HS and Founder VCFs by chromosome. 
    2. Run distinctive SNP filtering per chromosome. 
    3. Merge all chromosome-level results into one final VCF. 

Each step reports runtime and memory usage. 
"""

import os 
import argparse 
import subprocess 
import time 
import psutil 
from pathlib import Path

def run_cmd(cmd, label=""): 
    ''' Run a shell command, monitor runtime + memory. ''' 
    print(f"[INFO] Running: {cmd}") 
    start = time.time() 
    process = psutil.Process(os.getpid()) 

    # Launch command
    proc = subprocess.Popen(cmd, shell=True) 
    peak_mem = 0 
    while proc.poll() is None: 
        try: 
            mem = process.memory_info().rss / (1024 ** 2) # MB 
            peak_mem = max(peak_mem, mem) 
        except psutil.Error: 
            pass 
        time.sleep(0.1) 

    runtime = time.time() - start 
    if proc.returncode != 0: 
        raise RuntimeError(f"[ERROR] Command failed: {cmd}") 
    print(f"[INFO] {label} finished in {runtime:.2f}s | Peak memory: {peak_mem:.2f} MB\n") 

def split_vcfs(hs_vcf, founder_vcf, outdir, chroms, hs_prefix="chr", founder_prefix=""): 
    '''Split VCFs by chromosome.'''
    split_outputs = [] 
    for chrom in chroms: 
        hs_out = f"{outdir}/hs.chr{chrom}.vcf.gz" 
        founder_out = f"{outdir}/founder.chr{chrom}.vcf.gz" 
        run_cmd(f"bcftools view -r {hs_prefix}{chrom} {hs_vcf} -Oz -o {hs_out}", label=f"Split HS chr{chrom}") 
        run_cmd(f"bcftools view -r {founder_prefix}{chrom} {founder_vcf} -Oz -o {founder_out}", label=f"Split Founder chr{chrom}")
        run_cmd(f"bcftools index -f {hs_out}", label=f"Index HS chr{chrom}") 
        run_cmd(f"bcftools index -f {founder_out}", label=f"Index Founder chr{chrom}") 
        split_outputs.append((hs_out, founder_out)) 
    return split_outputs 

def process_chromosomes(split_outputs, outdir): 
    '''Run distinctive_snps.py per chromosome.''' 
    filtered_outputs = [] 
    for hs_vcf, founder_vcf in split_outputs: 
        chrom = Path(hs_vcf).stem.split('chr')[-1].replace('.vcf', '') 
        out_vcf = f"{outdir}/filtered_hs_chr{chrom}.vcf.gz" 
        cmd = ( 
               f"python3 distinctive_snps.py"
               f" --hs-vcf {hs_vcf}"
               f" --founder-vcf {founder_vcf}" 
               f" --chrom {chrom}" 
               f" --output-dir {outdir}" 
               f" --skip-snp-filtering" 
         )
        run_cmd(cmd, label=f"Process chr{chrom}") 
        filtered_outputs.append(out_vcf) 
    return filtered_outputs 

def merge_results(filtered_outputs, outdir): 
    '''Merge all chromosome-specific filtered results.'''
    final_vcf=f"{outdir}/filtered_all.vcf.gz" 
    inputs = " ".join(filtered_outputs) 
    run_cmd(f"bcftools concat -Oz -o {final_vcf} {inputs}", label=f"Concatenate all chromosomes") 
    run_cmd(f"bcftools index -f {final_vcf}", label=f"Index merged VCF") 
    return final_vcf 

# def merge_results(filtered_outputs, output_dir):
#     """
#     Merge per-chromosome filtered VCFs into a single file.
#     If sample sets are identical across files → bcftools concat
#     Otherwise → bcftools merge
#     """
#     final_vcf = os.path.join(output_dir, "filtered_all.vcf.gz")
#     inputs = " ".join(filtered_outputs)

#     # --- check sample sets ---
#     sample_sets = []
#     for vcf in filtered_outputs:
#         samples = subprocess.check_output(
#             f"bcftools query -l {vcf}", shell=True, text=True
#         ).strip().split("\n")
#         sample_sets.append(set(samples))

#     # If all sets are identical → concat
#     if all(s == sample_sets[0] for s in sample_sets):
#         print("[INFO] Using bcftools concat (sample sets identical across chromosomes)")
#         run_cmd(f"bcftools concat -Oz -o {final_vcf} {inputs}", 
#                 label="Concatenate all chromosomes")
#     else:
#         print("[INFO] Using bcftools merge (sample sets differ across chromosomes)")
#         run_cmd(f"bcftools merge -Oz -o {final_vcf} {inputs}", 
#                 label="Merge all chromosomes")

#     # Index the final file
#     run_cmd(f"bcftools index -f {final_vcf}", label="Index final VCF")
#     return final_vcf


def main(): 
    parser=argparse.ArgumentParser(description="Run SNP filtering pipeline with runtime/memory monitoring.") 
    parser.add_argument("--hs-vcf", required=True, help="Input HS VCF file (gzipped).") 
    parser.add_argument("--founder-vcf", required=True, help="Input Founder VCF file (gzipped).") 
    parser.add_argument("--chromosomes", nargs="+", required=True, help="List of chromosomes (e.g 1 2 3).") 
    parser.add_argument("--output-dir", required=True, help="Output directory.") 
    parser.add_argument("--hs-prefix", default="chr", help="Prefix for HS chromosomes (default: 'chr').") 
    parser.add_argument("--founder-prefix", default="", help="Prefix for Founder chromosomes (default:'').")
    args = parser.parse_args() 

    Path(args.output_dir).mkdir(parents=True, exist_ok=True) 

    # Step 1: Split 
    split_outputs = split_vcfs(args.hs_vcf, args.founder_vcf, args.output_dir, args.chromosomes, args.hs_prefix, args.founder_prefix) 

    # Step 2: Process 
    filtered_outputs = process_chromosomes(split_outputs, args.output_dir) 

    # Step 3: Merge 
    final_vcf = merge_results(filtered_outputs, args.output_dir)

    print(f"[INFO] Final merged VCF saved at: {final_vcf}") 

if __name__ == "__main__":
    main() 
