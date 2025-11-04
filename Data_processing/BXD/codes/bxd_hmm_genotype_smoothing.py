#!/usr/bin/env python3
"""
bxd_hmm_smoothing.py

Usage:
    python bxd_hmm_smoothing.py --input /path/to/BXD_no_meta.geno --out /path/to/BXD_hap_blocks.csv

What it does:
- Reads a genotype matrix with header columns:
    [chr_col, rs_col, cM_col, Mb_col, sample1, sample2, ...]
  where sample calls are 'B' (maternal), 'D' (paternal), 'H' (heterozygote), 'U' (unknown).
- For each chromosome and each sample, runs a 2-state HMM (states: B, D) using Viterbi
  with transition probabilities derived from genetic distance (cM).
- Collapses the Viterbi path into contiguous haplotype blocks and records start/end markers,
  physical (Mb) and genetic (cM) coordinates, number of markers, and a simple match fraction.
- Writes a CSV of haplotype blocks and prints a brief QC summary.

Dependencies:
    Python 3.8+
    pandas
    numpy
    (optional) tqdm for a progress bar

Install dependencies:
    pip install pandas numpy tqdm

Author: Felix Lisso
Date: 03-11-2025 
"""

import argparse
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

try:
    from tqdm import tqdm
except Exception:
    tqdm = lambda x: x  # no-progress fallback


# -------------------------
# Utilities / HMM settings
# -------------------------
# Emission probabilities: P(observed_call | hidden_state)
# These reflect a high-confidence match for the true founder, small probability for H or the opposite allele.
DEFAULT_EM = {
    'B': {'B': 0.99, 'H': 0.005, 'D': 0.004, 'U': 0.001},
    'D': {'D': 0.99, 'H': 0.005, 'B': 0.004, 'U': 0.001},
}

NEGLOG = lambda p: -math.log(max(p, 1e-12))


def recombination_fraction_from_cM(d_cM, min_r=1e-6, cap=0.5):
    """
    Approximate recombination fraction from genetic distance d in cM.
    For small distances, r ≈ d/100. We clip to (min_r, cap).
    """
    if d_cM is None or np.isnan(d_cM):
        return min_r
    r = abs(float(d_cM)) / 100.0
    if r <= 0:
        r = min_r
    return min(max(r, min_r), cap)


def clean_call(x):
    """Map raw cell to canonical {B, D, H, U}."""
    if pd.isna(x):
        return "U"
    s = str(x).strip().upper()
    if s.startswith("-"):
        s = s[1:]
    if s in {"B", "D", "H", "U"}:
        return s
    if len(s) > 0 and s[0] in {"B", "D", "H", "U"}:
        return s[0]
    return "U"


def viterbi(obs_list, cm_list, EM=DEFAULT_EM):
    """
    Viterbi algorithm for 2-state HMM (states 'B' and 'D').
    Inputs:
       obs_list: list of observations e.g. ['B','U','H','D',...]
       cm_list: numeric cM positions (same length). If missing, pass zeros or approximations.
    Returns:
       state_path: list of states ['B','B','D',...]
    """
    n = len(obs_list)
    if n == 0:
        return []

    # initial negative log-probabilities (assume equal prior 0.5)
    V_prev = {'B': NEGLOG(0.5) + NEGLOG(EM['B'].get(obs_list[0], 1e-12)),
              'D': NEGLOG(0.5) + NEGLOG(EM['D'].get(obs_list[0], 1e-12))}
    path = {'B': ['B'], 'D': ['D']}

    for i in range(1, n):
        d = cm_list[i] - cm_list[i - 1] if (cm_list[i] is not None and cm_list[i - 1] is not None) else 0.0
        r = recombination_fraction_from_cM(d)
        # transitions
        tBB, tBD = 1 - r, r
        tDD, tDB = 1 - r, r

        obs = obs_list[i]
        # compute negative log scores for each current state as min(prev + trans + emission)
        # For 'B' state at position i:
        score_via_B = V_prev['B'] + NEGLOG(tBB) + NEGLOG(EM['B'].get(obs, 1e-12))
        score_via_D = V_prev['D'] + NEGLOG(tDB) + NEGLOG(EM['B'].get(obs, 1e-12))
        if score_via_B <= score_via_D:
            best_B, best_prev_B = score_via_B, 'B'
        else:
            best_B, best_prev_B = score_via_D, 'D'

        # For 'D' state at position i:
        score_via_D = V_prev['D'] + NEGLOG(tDD) + NEGLOG(EM['D'].get(obs, 1e-12))
        score_via_B2 = V_prev['B'] + NEGLOG(tBD) + NEGLOG(EM['D'].get(obs, 1e-12))
        if score_via_D <= score_via_B2:
            best_D, best_prev_D = score_via_D, 'D'
        else:
            best_D, best_prev_D = score_via_B2, 'B'

        V_curr = {'B': best_B, 'D': best_D}
        path = {s: (path[best_prev_B] + ['B']) if s == 'B' else (path[best_prev_D] + ['D']) for s in ['B', 'D']}
        V_prev = V_curr

    final_state = 'B' if V_prev['B'] <= V_prev['D'] else 'D'
    return path[final_state]


# -------------------------
# Main processing
# -------------------------
def process_file(input_path, output_path, chr_col=None, rs_col=None, cm_col=None, mb_col=None, em=DEFAULT_EM, progress=True):
    # 1) Read file (auto-detect separator)
    print("Reading input:", input_path)
    df = pd.read_csv(input_path, sep=None, engine='python', dtype=str)
    print("Raw columns:", list(df.columns)[:10])

    # If user did not provide mapping, assume typical layout:
    # [chr, marker_id (rs), cM, Mb, sample...]
    if chr_col is None:
        chr_col = df.columns[0]
    if rs_col is None:
        rs_col = df.columns[1]
    if cm_col is None:
        cm_col = df.columns[2]
    if mb_col is None:
        mb_col = df.columns[3]

    sample_cols = list(df.columns[4:])
    print(f"Using mapping -> chr: {chr_col}, rs: {rs_col}, cM: {cm_col}, Mb: {mb_col}")
    print(f"Detected {len(sample_cols)} sample columns (first 6): {sample_cols[:6]}")

    # 2) Clean and prepare numeric columns
    df[cm_col] = pd.to_numeric(df[cm_col], errors='coerce')
    df[mb_col] = pd.to_numeric(df[mb_col], errors='coerce')

    # 3) Clean genotype calls
    for c in sample_cols:
        df[c] = df[c].apply(clean_call)

    # 4) Sort by chromosome and genetic map
    df_sorted = df.sort_values(by=[chr_col, cm_col], na_position='last').reset_index(drop=True)
    chromosomes = df_sorted[chr_col].unique().tolist()
    print("Chromosomes found:", chromosomes)

    # 5) Iterate chromosomes and samples, run Viterbi, collapse blocks
    blocks = []
    chrom_iter = tqdm(chromosomes) if progress else chromosomes
    for chrom in chrom_iter:
        chrom_df = df_sorted[df_sorted[chr_col] == chrom].reset_index(drop=True)
        if chrom_df.shape[0] == 0:
            continue
        n_markers = chrom_df.shape[0]
        marker_ids = chrom_df[rs_col].tolist()

        # Prepare cM list: forward/backfill then fallback to Mb or evenly spaced index-based cM
        cm_series = chrom_df[cm_col].fillna(method='ffill').fillna(method='bfill')
        if cm_series.isna().any():
            cm_series = cm_series.fillna(chrom_df[mb_col]).fillna(pd.Series(np.arange(n_markers) * 0.001))
        cm_list = cm_series.tolist()
        mb_list = chrom_df[mb_col].tolist()

        # For each sample
        for sample in chrom_df.columns[4:]:
            obs = chrom_df[sample].tolist()
            state_path = viterbi(obs, cm_list, EM=em)
            # collapse contiguous identical states to produce block entries
            start_idx = 0
            for i in range(1, n_markers + 1):
                if i == n_markers or state_path[i] != state_path[start_idx]:
                    end_idx = i - 1
                    assigned_state = state_path[start_idx]
                    block_obs = obs[start_idx:end_idx + 1]
                    match_frac = sum(1 for x in block_obs if x == assigned_state) / max(1, len(block_obs))
                    blocks.append({
                        'chr': chrom,
                        'sample': sample,
                        'start_marker': marker_ids[start_idx],
                        'end_marker': marker_ids[end_idx],
                        'start_mb': mb_list[start_idx],
                        'end_mb': mb_list[end_idx],
                        'start_cm': cm_list[start_idx],
                        'end_cm': cm_list[end_idx],
                        'num_markers': end_idx - start_idx + 1,
                        'state': assigned_state,
                        'match_fraction': match_frac
                    })
                    start_idx = i

    # 6) Save blocks table
    blocks_df = pd.DataFrame(blocks)
    out_cols = ['chr', 'sample', 'start_marker', 'end_marker',
                'start_mb', 'end_mb', 'start_cm', 'end_cm',
                'num_markers', 'state', 'match_fraction']
    blocks_df = blocks_df[out_cols]
    blocks_df.to_csv(output_path, index=False)
    print("Saved haplotype blocks to:", output_path)

    # 7) Quick QC prints
    print("Total blocks:", len(blocks_df))
    # Blocks per sample (top 10)
    summary = blocks_df.groupby('sample').size().reset_index(name='n_blocks').sort_values('n_blocks', ascending=False)
    print("Top 10 samples by number of blocks:")
    print(summary.head(10).to_string(index=False))

    return blocks_df


# -------------------------
# CLI
# -------------------------
def main():
    parser = argparse.ArgumentParser(description="BXD 2-state HMM smoothing & haplotype block generation")
    parser.add_argument('--input', '-i', required=True, help='Input genotype file (delimited text).')
    parser.add_argument('--out', '-o', required=True, help='Output CSV for haplotype blocks.')
    parser.add_argument('--no-progress', action='store_true', help='Disable progress bar.')
    # Optional explicit column names if different order
    parser.add_argument('--chr-col', help='Chromosome column name (default: first column).')
    parser.add_argument('--rs-col', help='Marker id column name (default: second column).')
    parser.add_argument('--cm-col', help='cM column name (default: third column).')
    parser.add_argument('--mb-col', help='Mb column name (default: fourth column).')

    args = parser.parse_args()
    input_path = Path(args.input)
    output_path = Path(args.out)

    if not input_path.exists():
        print("Input file not found:", input_path)
        sys.exit(2)

    process_file(input_path, output_path,
                 chr_col=args.chr_col, rs_col=args.rs_col, cm_col=args.cm_col, mb_col=args.mb_col,
                 progress=not args.no_progress)


if __name__ == '__main__':
    main()

