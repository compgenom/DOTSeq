#!/usr/bin/env python3
# coding: utf-8

"""
@author      CS Lim
@create date 2025-08-18 17:23:00
@modify date 2025-08-19 07:32:55
@desc        Generate ORF-centric, flatten GTF files for DOTSeq
"""

import os
import pickle
import sys
import argparse
import pandas as pd
import pyranges as pr
import csv
from riboss.orfs import orf_finder
from riboss.utils import flatten_orfs, bed6_to_bed12, bed12_to_gtf


def main():
    parser = argparse.ArgumentParser(description="Generate flattened GTF for DOTSeq and featureCounts.")
    parser.add_argument("--gtf", required=True, help="Path to GENCODE GTF file")
    parser.add_argument("--orf_finder", help="Path to RIBOSS orf_finder pickle file")
    parser.add_argument("--transcripts", help="Path to transcript FASTA file")
    parser.add_argument("--output", required=True, help="Output file prefix")
    parser.add_argument("--start_codon", default="ATG", help="Start codon to use")
    parser.add_argument("--stop_codon", action="store_true", help="Require stop codon")

    args = parser.parse_args()
    
    output_dir = os.path.dirname(args.output)
    if output_dir:
        if not os.path.exists(output_dir):
            try:
                os.makedirs(output_dir)
            except Exception as e:
                raise RuntimeError(
                    f"Failed to create the output directory '{output_dir}': {e}. "
                    "Please check that the path is valid and that you have write permissions."
                )
        elif not os.access(output_dir, os.W_OK):
            raise PermissionError(
                f"The output directory '{output_dir}' exists but is not writable. "
                "Please choose a different location or adjust the directory permissions."
            )

    if args.gtf and args.transcripts and not args.orf_finder:
        cds_range, df = orf_finder(args.gtf, args.transcripts, ncrna=False, outdir=None,
                                   start_codon=args.start_codon, stop_codon=args.stop_codon)
    elif args.orf_finder:
            df = pd.read_pickle(args.orf_finder)
    else:
        print("Error: --transcripts is required when using --gtf.", file=sys.stderr)
        sys.exit(1)

    orf_df, gt_map = flatten_orfs(df, args.gtf)
    
    # Drop overlapping genes with opposite strandness
    morfs = orf_df[orf_df.Score == "mORF"].drop_duplicates("Name", keep = False)
    orf_df = pd.merge(orf_df, morfs[["Name"]])

    bed12_df = bed6_to_bed12(orf_df)
    dotseq_gtf = bed12_to_gtf(bed12_df, gt_map)

    dotseq_gtf.to_csv(f"{args.output}.gtf", sep="\t", header=False, index=False,
                      quoting=csv.QUOTE_NONE, quotechar='')
    orf_df.to_csv(f"{args.output}.bed", sep="\t", header=False, index=False)
    
if __name__ == "__main__":
    main()
