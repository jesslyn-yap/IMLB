import pyranges as pr
import os
import pandas as pd
import numpy as np

def get_genome_sizes(fai_path):
    if not os.path.exists(fai_path):
        raise FileNotFoundError(f"FAI index file not found: {fai_path}\n")
    
    genome_df = pd.read_csv(
        fai_path,
        sep="\t",
        header=None,
        usecols=[0, 1],
        names=["Chromosome", "Length"]
    )
    genome_df["Length"] = genome_df["Length"].astype(int)

    # Keep only chr1–22
    genome_df = genome_df[genome_df["Chromosome"].str.match(r"^chr([1-9]|1[0-9]|2[0-2])$")]

    genome_sizes = genome_df.set_index("Chromosome")["Length"].to_dict()

    return genome_sizes

def get_bins(fai_path, out_folder):
    genome_sizes = get_genome_sizes(fai_path)

    bins = pd.DataFrame({
        "Chromosome": [],
        "Start": [],
        "End": []
    })

    for chrom, length in genome_sizes.items():
        starts = list(range(0, length, 200))
        ends = [s+200 for s in starts if s+200 <= length]
        starts = starts[:len(ends)]

        chrom_df = pd.DataFrame({
            "Chromosome": chrom,
            "Start": starts,
            "End": ends
        })

        bins = pd.concat([bins, chrom_df], ignore_index = True)

    os.makedirs(out_folder, exist_ok=True)
    out_path = os.path.join(out_folder, "genome_bins.bed")
    bins.to_csv(out_path, sep='\t', header=False, index=False)

    print("Bins generated")
    return bins

def create_labeled_data(tf, fai_path, out_folder):
    """
    Generate tsv file containing labels for genomic bins:
    - Positive: bins overlapping the combined ChIP-seq peaks.
    - Ambiguous: bins overlapping at least one replicate but not in the combined peaks.
    - Negative: all other bins (bins that do not overlap any peaks).
    Ambiguous bins are discarded from the final dataset.

    Args:
        tf: Interested transcription factor
        fai_path: absolute path to .fai file
        out_folder: absolute path to output folder
    """
    bins_file = os.path.join(out_folder, "genome_bins.bed")

    if not os.path.exists(bins_file):
        bins_df = get_bins(fai_path, out_folder)
        bins = pr.PyRanges(bins_df)
    else:
        bins = pr.read_bed(bins_file)

    label_file = os.path.join(out_folder, f"{tf}_labelled.bed")

    if not os.path.exists(label_file):
        tf_combined = pr.read_bed(f"{tf}-bed/{tf}-bed_combined.bed")
        tf_1 = pr.read_bed(f"{tf}-bed/{tf}-bed_1.bed")
        tf_2 = pr.read_bed(f"{tf}-bed/{tf}-bed_2.bed")

        pos = bins.overlap(tf_combined)
        pos_df = pos.df.copy()
        pos_df["label"] = 1

        neg = bins.overlap(tf_combined, invert=True)
        neg = neg.overlap(tf_1, invert=True)
        neg = neg.overlap(tf_2, invert=True)
        neg_df = neg.df.copy()
        neg_df["label"] = 0

        all_df = pd.concat([pos_df, neg_df], ignore_index=True)
        all_df = all_df.sample(frac=1, random_state=42).reset_index(drop=True)

        all_df.to_csv(f"{out_folder}/{tf}_labels.tsv", sep="\t", index=False)
    else:
        print(f"{tf} labeled file already exists, skipping...")

    print("Labels generated")

