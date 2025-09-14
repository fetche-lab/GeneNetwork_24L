#!/usr/bin/env python3

import pandas as pd
import numpy as np
import os
import subprocess 
import gzip
import csv
import argparse
import time
import psutil


def normalize_chrom(chrom):
    return str(chrom).lower().replace("chr", "")


def get_vcf_positions(vcf_file, chrom="1", min_call_rate=0.95):
    positions = set()
    with gzip.open(vcf_file, 'rt') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            if row[0].startswith('#'):
                continue
            norm_chrom = normalize_chrom(row[0])
            if norm_chrom != chrom:
                continue
            pos = int(row[1])
            format_field = row[8]
            if 'GT' not in format_field:
                continue
            gt_index = format_field.split(':').index('GT')
            gts = [sample.split(':')[gt_index] for sample in row[9:]]
            call_count = sum(1 for gt in gts if gt not in ['./.', '.|.'])
            call_rate = call_count / len(gts) if len(gts) > 0 else 0
            if call_rate >= min_call_rate:
                positions.add(pos)
    return positions


def load_genotypes_at_positions(vcf_file, positions, is_founder=False):
    genotypes = {pos: [] if is_founder else {} for pos in positions}
    with gzip.open(vcf_file, 'rt') as f:
        reader = csv.reader(f, delimiter='\t')
        header = next(reader)
        while header[0].startswith('#'):
            header = next(reader)
        samples = header[9:]
        for row in reader:
            if row[0].startswith('#'):
                continue
            pos = int(row[1])
            if pos not in positions:
                continue
            ref = row[3]
            alt = row[4].split(',')
            format_field = row[8]
            gt_index = format_field.split(':').index('GT')
            for i, sample_field in enumerate(row[9:]):
                gt_str = sample_field.split(':')[gt_index]
                if gt_str == './.' or gt_str == '.|.':
                    allele = np.nan
                else:
                    gt = gt_str.split('|') if '|' in gt_str else gt_str.split('/')
                    if len(gt) != 2:
                        allele = np.nan
                    else:
                        mat_gt, pat_gt = int(gt[0]), int(gt[1])
                        mat_allele = ref if mat_gt == 0 else alt[mat_gt - 1]
                        pat_allele = ref if pat_gt == 0 else alt[pat_gt - 1]
                        if is_founder:
                            allele = mat_allele  # homozygous founders
                        else:
                            allele = f"{mat_allele}|{pat_allele}"
                if is_founder:
                    genotypes[pos].append(allele)
                else:
                    genotypes[pos][samples[i]] = allele
    return genotypes


def load_and_intersect_vcfs(hs_vcf_file, founder_vcf_file, chrom="1", min_call_rate=0.95,
                             snp_call_rate=1.0, skip_snp_filtering=False):
    hs_positions = get_vcf_positions(hs_vcf_file, chrom, min_call_rate)
    founder_positions = get_vcf_positions(founder_vcf_file, chrom, min_call_rate)
    common_positions = sorted(hs_positions & founder_positions)
    print(f"[INFO] Intersecting SNP positions: {len(common_positions)}")

    hs_genotypes = load_genotypes_at_positions(hs_vcf_file, common_positions, is_founder=False)
    founder_genotypes = load_genotypes_at_positions(founder_vcf_file, common_positions, is_founder=True)

    snp_info = pd.DataFrame({'chr': chrom, 'bp_pos': common_positions, 'pos_mb': [p / 1e6 for p in common_positions]})
    hs_df = pd.DataFrame([hs_genotypes[pos] for pos in common_positions])
    founder_df = pd.DataFrame([founder_genotypes[pos] for pos in common_positions])

    if not skip_snp_filtering and snp_call_rate > 0:
        hs_call_rate = hs_df.notna().mean(axis=1)
        founder_call_rate = founder_df.notna().mean(axis=1)
        keep_snps = (hs_call_rate >= snp_call_rate) & (founder_call_rate >= snp_call_rate)

        dropped_snps = (~keep_snps).sum()
        print(f"[INFO] Dropped {dropped_snps} SNPs due to call rate filtering")

        snp_info = snp_info[keep_snps].reset_index(drop=True)
        hs_df = hs_df[keep_snps].reset_index(drop=True)
        founder_df = founder_df[keep_snps].reset_index(drop=True)
    else:
        print("[INFO] Skipping SNP-level call rate filtering")

    print(f"[INFO] SNPs after filtering: {len(snp_info)}")
    return hs_df, founder_df, snp_info


def identify_distinctive_positions(hs_genotypes, founder_genotypes, snp_info, min_diff_pct=80):
    n_snps, n_samples = hs_genotypes.shape
    n_founders = founder_genotypes.shape[1]
    min_diff_founders = int(n_founders * (min_diff_pct / 100))
    distinctive_pos = set()

    for start in range(0, n_snps, 10000):
        end = min(start + 10000, n_snps)
        hs_chunk = hs_genotypes.iloc[start:end].values
        founder_chunk = founder_genotypes.iloc[start:end].values.astype(str)

        for sample_idx in range(n_samples):
            sample_col = hs_chunk[:, sample_idx]
            valid_mask = ~pd.isna(sample_col)
            if not np.any(valid_mask):
                continue

            mat_alleles = np.array([x.split('|')[0] if pd.notna(x) else np.nan for x in sample_col])
            pat_alleles = np.array([x.split('|')[1] if pd.notna(x) else np.nan for x in sample_col])

            for alleles in [mat_alleles, pat_alleles]:
                valid_alleles = alleles[valid_mask]
                if len(valid_alleles) == 0:
                    continue
                diff_counts = np.sum(founder_chunk[valid_mask] != valid_alleles[:, np.newaxis], axis=1)
                distinctive_mask = diff_counts >= min_diff_founders
                if np.any(distinctive_mask):
                    local_indices = np.where(distinctive_mask)[0]
                    global_indices = start + np.where(valid_mask)[0][local_indices]
                    distinctive_pos.update(snp_info.iloc[global_indices]['bp_pos'])

    print(f"[INFO] Distinctive SNPs found: {len(distinctive_pos)}")
    return distinctive_pos


def create_filtered_vcf(hs_vcf_file, distinctive_pos, output_vcf_file):
    with gzip.open(hs_vcf_file, 'rt') as infile, gzip.open(output_vcf_file, 'wt') as outfile:
        for line in infile:
            if line.startswith('#'):
                outfile.write(line)
                continue
            fields = line.strip().split('\t')
            pos = int(fields[1])
            if pos in distinctive_pos:
                outfile.write(line)

    print(f"[INFO] Filtered VCF saved to {output_vcf_file}")


def main():
    parser = argparse.ArgumentParser(description="Identify distinctive SNPs in HS vs Founder VCFs")
    parser.add_argument("--hs-vcf", required=True, help="VCF file for HS samples (gzipped)")
    parser.add_argument("--founder-vcf", required=True, help="VCF file for Founder samples (gzipped)")
    parser.add_argument("--chrom", required=True, help="Chromosome to process")
    parser.add_argument("--output-dir", default="output", help="Directory for outputs")
    parser.add_argument("--min-diff-pct", type=float, default=80, help="Min percent founders differing to call distinctive SNP")
    parser.add_argument("--snp-call-rate", type=float, default=1.0, help="Min SNP call rate (0-1). Use 0 to disable")
    parser.add_argument("--skip-snp-filtering", action="store_true", help="Skip SNP-level call rate filtering")

    args = parser.parse_args()

    start_time = time.time()

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    # Load and intersect SNPs
    hs_positions = get_vcf_positions(args.hs_vcf, args.chrom)
    founder_positions = get_vcf_positions(args.founder_vcf, args.chrom)
    common_positions = sorted(hs_positions & founder_positions)

    if len(common_positions) == 0:
        print(f"[INFO] No overlapping SNPs found on chromosome {args.chrom}. Skipping...")

        out_vcf = os.path.join(args.output_dir, f"filtered_hs_chr{args.chrom}.vcf.gz")

        # Write an empty but valid VCF
        tmp_vcf = out_vcf.replace(".gz", "")
        sample_ids = [
    "56269-4", "56272-3", "56275-4", "57215-3", "57241-4", "57241-5",
    "57604-6", "57810-3", "57909-4", "57914-3", "58041-6", "58095-3",
    "58096-4", "58096-6", "58105-4", "58125-3", "58125-4", "58293-3",
    "58293-4", "58294-4", "58295-4", "58301-3", "58301-4", "58337-3",
    "58342-3", "58343-4", "58564-4", "58566-3", "58566-4", "58570-3",
    "58572-4", "58765-4", "58769-3", "58771-4", "58774-4", "58779-3",
    "58779-4", "59104-4", "59108-3", "59114-3", "59114-4", "60213-3",
    "60214-3", "60214-5", "60444-3", "60444-5", "60445-3", "60641-3",
    "60644-4", "60650-3", "61038-4", "61038-5", "61040-5", "61158-4",
    "61160-3", "61160-4", "61163-4", "61172-3", "61357-4", "61638-4",
    "61638-5", "61638-6", "61641-3", "61643-5", "61644-4", "61645-3",
    "61646-3", "61646-4", "61821-6", "62044-5", "62046-4", "62288-3",
    "62289-4", "62296-3", "62296-4", "62298-4", "62302-4", "62309-3",
    "62309-4", "62310-4", "62312-4", "62315-4", "62316-4", "62317-4",
    "62318-4", "62321-3", "62322-4", "62322-6", "62323-3", "62324-6",
    "62329-4", "62333-4", "62334-4", "62335-3", "62335-4", "62338-3",
    "62339-3", "62339-4", "62343-4", "62347-5", "62347-6", "62349-3",
    "62349-6", "62352-4", "62360-5", "62372-3", "62378-3", "62378-4",
    "62380-4", "62381-5", "62381-6", "62397-4", "62407-3", "62407-6",
    "62408-3", "62411-4", "62414-5", "62414-6", "62415-3", "62417-3",
    "62418-3", "62420-3", "62422-3", "62422-4", "62425-4", "62426-3",
    "62432-3", "62436-3", "62436-4", "62437-4", "62438-3", "62438-5",
    "62443-4", "62449-3", "62471-3", "62485-4", "62495-3", "62496-4",
    "62497-4", "62506-4", "62507-3", "62507-4", "62511-3", "62514-4",
    "62516-4", "62561-3", "62561-4", "62562-5", "62562-6", "62594-3",
    "62594-6", "62595-3", "62596-3", "62597-4", "62607-3", "62608-4",
    "62608-5", "62608-6", "62611-3", "62612-4", "62613-3", "62615-3",
    "62618-4", "62631-3", "62633-3", "62633-4", "62634-3", "62635-4",
    "62636-4", "62641-4", "62642-3", "62642-4", "62647-3", "62649-3",
    "62651-4", "62652-3", "62652-4", "62653-3", "62653-4", "62655-3",
    "62658-3", "62661-4", "62663-3", "62663-4", "62667-4", "62669-3",
    "62669-5", "62675-4", "62677-3", "62678-4", "62678-6", "62679-3",
    "62680-3", "62683-4", "62684-3", "62684-4", "62685-3", "62686-4",
    "62692-3", "62695-4", "62700-3", "62718-4", "62720-3", "62724-3",
    "62725-3", "62728-3", "62728-4", "62731-4", "62732-3", "62733-5",
    "62734-4", "62762-4", "62763-3", "62765-4", "62769-3", "62769-6",
    "62785-4", "62789-4", "62792-4", "62797-3", "62798-3", "62799-4",
    "62800-3", "62818-6", "62819-3", "62819-4", "62854-4", "62857-4",
    "62864-4", "62868-3", "62872-4", "62873-5", "62875-3", "62882-4",
    "62883-3", "62884-3", "62886-3", "62887-3", "62888-3", "62892-3",
    "62892-4", "62895-3", "62896-3", "62896-4", "62897-4", "62900-3",
    "62901-3", "62901-4", "62901-5", "62909-4", "62910-4", "62911-3",
    "62911-4", "62912-4", "62913-3", "62913-6", "62919-4", "62920-5",
    "62921-4", "62922-4", "62934-3", "62935-3", "62936-4", "62938-3",
    "62943-3", "62943-4", "62961-3", "62966-3", "62968-3", "62975-3",
    "62978-6", "62980-3", "62988-6", "62995-4", "62996-3", "62996-4",
    "63000-3", "63000-4", "63002-4", "63006-4", "63008-4", "63008-5",
    "63014-4", "63016-4", "63032-3", "63032-4", "63032-5", "63056-3",
    "63906-4", "63910-4", "63912-3", "63912-4", "63913-4", "63939-3",
    "63939-4", "63940-3", "63940-4", "63941-3", "63941-4", "63942-3",
    "63942-4", "63942-5", "63944-3", "63944-5", "63947-3", "63947-4",
    "63948-3", "63948-4", "63949-3", "63952-3", "63952-4", "63953-3",
    "63956-4", "63956-5", "63959-3", "63959-4", "63961-3", "63961-4",
    "63966-3", "63966-4", "63968-4", "63972-3", "63972-5", "63974-4",
    "63984-3", "63984-4", "63990-4", "63991-3", "63997-3", "64000-3",
    "64000-4", "64001-3", "64003-3", "64003-4", "64003-5", "64009-3",
    "64011-3", "64016-4", "64016-6", "64023-4", "64024-3", "64024-4",
    "64027-3", "64027-4", "64034-3", "64038-3", "64038-4", "64039-3",
    "64043-3", "64045-3", "64045-4", "64057-4", "64060-3", "64070-4",
    "64104-3", "64105-3", "64105-4", "64106-3", "64106-5", "64108-3",
    "64109-3", "64109-4", "64110-3", "64111-3", "64112-3", "64113-3",
    "64113-4", "64116-4", "64119-3", "64131-4", "64134-3", "64134-5",
    "64134-6", "64142-4", "64146-3", "64146-4", "64148-3", "64148-4",
    "64149-3", "64152-4", "64154-4", "64159-3", "64159-4", "64160-2",
    "64162-3", "64162-4", "64164-3", "64166-4", "64167-3", "64168-4",
    "64170-3", "64170-4", "64171-3", "64177-4", "64180-3", "64183-3",
    "64185-4", "64186-4", "64194-3", "64199-4", "64199-5", "64205-3",
    "64213-3", "64213-4", "64226-4", "64227-3", "64230-4", "64230-5",
    "64231-3", "64231-4", "64235-3", "64235-5", "64238-3", "64238-4",
    "64240-3", "64240-4", "64244-4", "64246-4", "64246-5", "64247-3",
    "64248-3", "64249-4", "64251-3", "64257-4", "64262-4", "64267-3",
    "64293-3", "64294-4", "64295-4", "64297-3", "64297-5", "64297-6",
    "64299-4", "64304-1", "64305-5", "64305-6"
]

        with open(tmp_vcf, "w") as f:
            f.write("##fileformat=VCFv4.2\n")
            header_cols =["#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"] + sample_ids
            f.write("\t".join(header_cols) + "\n")

        subprocess.run(f"bgzip -f {tmp_vcf}", shell=True, check=True)
        subprocess.run(f"bcftools index -f {out_vcf}", shell=True, check=True)

        return  # exit gracefully so downstream doesn’t crash

    # If we reach here, continue with filtering + distinctive SNPs
    hs_df, founder_df, snp_info = load_and_intersect_vcfs(
        args.hs_vcf, args.founder_vcf, args.chrom,
        snp_call_rate=args.snp_call_rate,
        skip_snp_filtering=args.skip_snp_filtering
    )

    distinctive_pos = identify_distinctive_positions(hs_df, founder_df, snp_info, min_diff_pct=args.min_diff_pct)
    output_vcf_file = os.path.join(args.output_dir, f"filtered_hs_chr{args.chrom}.vcf.gz")
    create_filtered_vcf(args.hs_vcf, distinctive_pos, output_vcf_file)

    end_time = time.time()
    elapsed = end_time - start_time
    mem_used = psutil.Process(os.getpid()).memory_info().rss / (1024**2)
    print(f"[INFO] Runtime: {elapsed:.2f} sec, Memory used: {mem_used:.2f} MB")

    # if len(common_positions) == 0:
    #     print(f"[INFO] No overlapping SNPs found on chromosome {chrom}. Skipping...")

    #     out_vcf = os.path.join(output_dir, f"filtered_hs_chr{chrom}.vcf.gz")

    #     # Write an empty but valid VCF
    #     tmp_vcf = out_vcf.replace(".gz", "")
    #     with open(tmp_vcf, "w") as f:
    #       f.write("##fileformat=VCFv4.2\n")
    #       f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

    #    subprocess.run(f"bgzip -f {tmp_vcf}", shell=True, check=True)
    #    subprocess.run(f"bcftools index -f {out_vcf}", shell=True, check=True)

    #    return  # exit gracefully so downstream doesn’t crash


    # if not os.path.exists(args.output_dir):
    #     os.makedirs(args.output_dir)

    # hs_df, founder_df, snp_info = load_and_intersect_vcfs(
    #     args.hs_vcf, args.founder_vcf, args.chrom,
    #     snp_call_rate=args.snp_call_rate,
    #     skip_snp_filtering=args.skip_snp_filtering
    # )

    # distinctive_pos = identify_distinctive_positions(hs_df, founder_df, snp_info, min_diff_pct=args.min_diff_pct)
    # output_vcf_file = os.path.join(args.output_dir, f"filtered_hs_chr{args.chrom}.vcf.gz")
    # create_filtered_vcf(args.hs_vcf, distinctive_pos, output_vcf_file)

    # end_time = time.time()
    # elapsed = end_time - start_time
    # mem_used = psutil.Process(os.getpid()).memory_info().rss / (1024**2)
    # print(f"[INFO] Runtime: {elapsed:.2f} sec, Memory used: {mem_used:.2f} MB")


if __name__ == "__main__":
    main()

