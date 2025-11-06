#!/usr/bin/env python3
"""
bxd_select_representative_markers.py

Usage:
    python bxd_select_representative_markers.py \
        --geno BXD_no_meta.geno \
        --blocks BXD_hap_blocks.csv \
        --out-geno BXD_reduced.geno \
        --out-markers selected_markers.csv

Goal:
    From the original genotype matrix and smoothed haplotype blocks,
    pick the most representative marker(s) per block and create a reduced
    version of the genotype file.

Dependencies:
    pandas, numpy
"""

import argparse
import pandas as pd
import numpy as np

def clean_call(x):
    if pd.isna(x): return "U"
    s = str(x).strip().upper()
    if s.startswith("-"): s = s[1:]
    return s if s in {"B","D","H","U"} else "U"

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--geno", required=True)
    parser.add_argument("--blocks", required=True)
    parser.add_argument("--out-geno", required=True)
    parser.add_argument("--out-markers", required=True)
    args = parser.parse_args()

    print(f"Reading genotype file: {args.geno}")
    geno = pd.read_csv(args.geno, sep=None, engine="python", dtype=str)
    chr_col, rs_col, cm_col, mb_col = geno.columns[:4]
    sample_cols = geno.columns[4:]
    print(f"Loaded {geno.shape[0]} markers × {len(sample_cols)} samples")

    geno[sample_cols] = geno[sample_cols].applymap(clean_call)
    geno[mb_col] = pd.to_numeric(geno[mb_col], errors="coerce")
    geno[cm_col] = pd.to_numeric(geno[cm_col], errors="coerce")

    print(f"Reading haplotype blocks: {args.blocks}")
    blocks = pd.read_csv(args.blocks)

    selected = []
    for idx, blk in blocks.iterrows():
        chrom = str(blk["chr"])
        s_mb, e_mb = blk["start_mb"], blk["end_mb"]
        state = blk["state"]
        block_df = geno[(geno[chr_col].astype(str)==chrom) &
                        (geno[mb_col]>=s_mb) & (geno[mb_col]<=e_mb)]
        if block_df.empty:
            continue

        # Compute marker-level metrics
        calls = block_df[sample_cols]
        callrate = (calls != "U").sum(axis=1) / len(sample_cols)
        concordance = (calls.apply(lambda r: (r==state).sum(), axis=1)) / len(sample_cols)
        fracB = (calls == "B").sum(axis=1) / len(sample_cols)
        fracD = (calls == "D").sum(axis=1) / len(sample_cols)
        informative = 2*np.minimum(fracB, fracD)   # 0–1
        score = 0.45*concordance + 0.30*callrate + 0.20*informative
        block_df = block_df.assign(callrate=callrate,
                                   concordance=concordance,
                                   informative=informative,
                                   score=score)
        best = block_df.loc[score.idxmax()]
        selected.append(best)

    sel_df = pd.DataFrame(selected)
    print(f"Selected {len(sel_df)} representative markers out of {geno.shape[0]}")

    # Save selected marker list
    marker_cols = [chr_col, rs_col, cm_col, mb_col, "callrate", "concordance", "informative", "score"]
    sel_df[marker_cols].to_csv(args.out_markers, index=False)

    # Write reduced genotype file (same format as input)
    reduced = geno[geno[rs_col].isin(sel_df[rs_col])]
    reduced.to_csv(args.out_geno, index=False)
    print(f"Wrote reduced genotype file: {args.out_geno}")
    print("Done.")

if __name__ == "__main__":
    main()

