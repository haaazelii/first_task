# This program analyzes protein alignment results from ssearch36 and domain annotations
# from CDD to identify conserved regions across proteins. It highlights regions shared
# by a user-defined fraction of proteins.

import pandas as pd
import re
import argparse
import matplotlib.pyplot as plt
import os

# sample command: python3 first_task.py -ssi () -cddi () -cdde () -cov () -outdir ()
parser = argparse.ArgumentParser()
parser.add_argument("-ssi", required=True, help="Path to ssearch input file.")
parser.add_argument("-cddi", required=True, help="Path to CDD input file.")
parser.add_argument("-cdde", required=True, type=float, help="E‑value cutoff for CDD hits (default 1e‑5)")
parser.add_argument("-cov", required=True, type=float, help="Coverage cutoff for shaded region")
parser.add_argument("-outdir", required=True, help="Directory to save PNG graphs")
args = parser.parse_args()

# for ssearch input to output
ssearch_input = args.ssi
query_lengths = {}
filtered_lines = []

with open(ssearch_input, "r") as file:
    for line in file:
        # store query length
        match = re.match(r"# Query: (\S+)  - (\d+) aa", line)
        if match:
            query_id = match.group(1)
            query_lengths[query_id] = int(match.group(2))
        # store only non-header lines
        elif not line.startswith("#"):
            filtered_lines.append(line)

columns = [
    "Query ID", "Subject ID", "% Identity", "Alignment Length", "Mismatches", "Gap Opens", 
    "Left Pos", "Right Pos", "Subject Start", "Subject End", "E-value", "Bit Score"
]

data = [line.split("\t") for line in filtered_lines]
df = pd.DataFrame(data, columns=columns)
count = len(pd.concat([df["Query ID"], df["Subject ID"]]).drop_duplicates())

# Convert numeric columns
num_cols = [
    "% Identity", "Alignment Length", "Mismatches", "Gap Opens", 
    "Left Pos", "Right Pos", "Subject Start", "Subject End", "E-value", "Bit Score"
]
for col in num_cols:
    df[col] = pd.to_numeric(df[col])

# Remove self comparisons
df_filtered = df[df["Query ID"] != df["Subject ID"]].copy()

# Compute Qcov
df_filtered["Query Coverage"] = df_filtered.apply(
    lambda row: (row["Right Pos"] - row["Left Pos"] + 1) / query_lengths.get(row["Query ID"]), axis=1
)

df_filtered["Query Length"] = df_filtered["Query ID"].map(query_lengths)
df_filtered["Subject Length"] = df_filtered["Subject ID"].map(query_lengths)

df_comparison = df_filtered[[
    "Query ID", "Subject ID", "Query Length", "Subject Length",
    "Left Pos", "Right Pos", "Subject Start", "Subject End", 
    "Query Coverage", "Bit Score", "E-value"
]]

df = df_comparison.iloc[:, :-3]

# for CDD input to output
cdd_cols = ["Query ID","query length","CDD ID","CDD length","evalue","bit score",
            "q. start","q. end","% query coverage per hsp","s. start","s. end"]

cdd_rows = []
with open(args.cddi) as fh:
    for ln in fh:
        if ln.startswith('#'):
            continue
        parts = ln.rstrip('\n').split('\t')
        if len(parts) >= 11:
            cdd_rows.append(parts[:11])

cdd_df = pd.DataFrame(cdd_rows, columns=cdd_cols)
for col in ["evalue","q. start","q. end"]:
    cdd_df[col] = pd.to_numeric(cdd_df[col], errors='coerce')

# e-value filter
cdd_df = cdd_df[cdd_df['evalue'] < args.cdde]
cdd_df = cdd_df.sort_values(by=["Query ID", "q. start"])

# start graphing
# extract the reversed match for graphing
reversed_df = df.copy()
reversed_df["Query ID"] = df["Subject ID"]
reversed_df["Subject ID"] = df["Query ID"]
reversed_df["Query Length"] = df["Subject Length"]
reversed_df["Subject Length"] = df["Query Length"]
reversed_df["Left Pos"] = df["Subject Start"]
reversed_df["Right Pos"] = df["Subject End"]
reversed_df["Subject Start"] = df["Left Pos"]
reversed_df["Subject End"] = df["Right Pos"]

# combine and sort
df_all = pd.concat([df, reversed_df], ignore_index=True)
df_all = df_all.drop_duplicates()
df_all = df_all.sort_values(by=["Query ID", "Left Pos"])

out_dir = args.outdir
os.makedirs(out_dir, exist_ok=True)

# for each query
for query in df_all["Query ID"].unique():
    subset = df_all[df_all["Query ID"] == query]
    length = subset["Query Length"].iloc[0]
    
    # Filter CDD annotations for this query
    cdd_subset = cdd_df[cdd_df["Query ID"] == query]

    fig, ax = plt.subplots(figsize=(10, 1.5 + 0.5 * (len(subset) + len(cdd_subset))))

    # Calculate coverage
    coverage = [set() for _ in range(length + 1)]
    for row in subset.itertuples():
        subj = row._2
        start = max(1, row._5)
        end = min(length, row._6)
        for i in range(start, end + 1):
            coverage[i].add(subj)

    # Identify significant shared regions using coverage
    sig_regions = []
    in_region = False
    region_start = None

    for i in range(1, length + 1):
        if (len((coverage[i])) / (count - 1))>= args.cov:
            if not in_region:
                region_start = i
                in_region = True
        else:
            if in_region:
                region_end = i - 1
                in_region = False

                # Extend region if it overlaps with a CDD and that CDD has wider boundaries
                for cdd in cdd_subset.itertuples():
                    cdd_start = cdd._7
                    cdd_end = cdd._8

                    # If region start is inside CDD and CDD starts earlier → extend start
                    if cdd_start <= region_start <= cdd_end and cdd_start < region_start:
                        region_start = cdd_start

                    # If region end is inside CDD and CDD ends later → extend end
                    if cdd_start <= region_end <= cdd_end and cdd_end > region_end:
                        region_end = cdd_end
                sig_regions.append((region_start, region_end))
    
    # Catch final region
    if in_region:
        region_end = length
        for cdd in cdd_subset.itertuples():
            cdd_start = cdd._7
            cdd_end = cdd._8
            if abs(region_start - cdd_end) <= 1:
                region_start = min(region_start, cdd_start)
            if abs(region_end - cdd_start) <= 1:
                region_end = max(region_end, cdd_end)
        sig_regions.append((region_start, region_end))

    # Plot shaded regions
    for start, end in sig_regions:
        ax.axvspan(start, end, color='red', alpha=0.2)
        ax.text((start + end)/2, len(subset) + len(cdd_subset) + 1.5, f"{start}-{end}", 
                ha='center', va='bottom', fontsize=8, color='red')

    # plot query sequence
    ax.hlines(
        y=len(subset) + len(cdd_subset) + 1,
        xmin=1,
        xmax=length,
        colors='black',
        linewidth=5,
    )
    ax.text(length + length*0.1, len(subset) + len(cdd_subset) + 1, query, va='center')

    # plot CDD hits
    for i, row in enumerate(reversed(list(cdd_subset.itertuples())), 1):
        y_pos = len(subset) + i
        ax.hlines(
            y=y_pos, 
            xmin=row._7, 
            xmax=row._8, 
            colors='green', 
            linewidth=5
        )
        ax.text(length + length*0.1, y_pos, row._3, va='center', color='black')

    # plot subject matches
    for i, row in enumerate(reversed(list(subset.itertuples())), 1):
        left_pos = row._5 
        right_pos = row._6  
        subject = row._2 
        ax.hlines(
            y=i, 
            xmin=left_pos, 
            xmax=right_pos, 
            colors='blue', 
            linewidth=5
        )
        ax.text(length + length*0.1, i, subject, va='center')

    ax.set_ylim(0, len(subset) + len(cdd_subset) + 2)
    ax.set_xlim(0, length + length*0.50)
    ax.set_xlabel("Query Sequence Position")
    ax.set_yticks([])
    ax.set_title(f"Query: {query} (coverage >= {args.cov}, n = {count})")

    fig.tight_layout()
    fig.savefig(f"{out_dir}/{query.replace('/', '_')}_alignment.png")
    plt.close()