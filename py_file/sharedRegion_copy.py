# This program analyzes protein alignment results from ssearch36 and domain annotations
# from CDD to identify conserved regions across proteins. It highlights regions shared
# by a user-defined fraction of proteins.

# for ssearch output and rescued cdd output

# count part of cdd not all

import pandas as pd
import re
import argparse
import matplotlib.pyplot as plt
import os

# sample command: python3 sharedRegion.py -ssi () -cddi () -cdde () -vote () -outdir ()
parser = argparse.ArgumentParser()
parser.add_argument("-ssi", required=True, help="Path to ssearch input file.")
parser.add_argument("-cddi", required=True, help="Path to CDD input file.")
parser.add_argument("-cdde", required=True, type=float, help="E‑value cutoff for CDD hits (default 1e‑5)")
parser.add_argument("-vote", required=True, type=float, help="vote cutoff for shaded region")
parser.add_argument("-outdir", required=True, help="Directory to save PNG graphs")
args = parser.parse_args()

# ssearch input
length_records = []
ssearch_input = args.ssi
filtered_lines = []

with open(ssearch_input, "r") as file:
    for line in file:
        match = re.match(r"# Query:\s+(\S+)\s+-\s+(\d+)\s+aa", line)
        if match:
            length_records.append((match.group(1), int(match.group(2))))
        elif not line.startswith("#"):
            filtered_lines.append(line)

length_df = pd.DataFrame(length_records, columns=["ID", "Length"])
length_lookup = length_df.set_index("ID")["Length"].to_dict()

columns = [
    "Query ID", "Subject ID", "% Identity", "Alignment Length", "Mismatches", "Gap Opens", 
    "Left Pos", "Right Pos", "Subject Start", "Subject End", "E-value", "Bit Score"
]

data = [line.split("\t") for line in filtered_lines]
df = pd.DataFrame(data, columns=columns)
count = len(pd.concat([df["Query ID"], df["Subject ID"]]).drop_duplicates())

num_cols = [
    "% Identity", "Alignment Length", "Mismatches", "Gap Opens", 
    "Left Pos", "Right Pos", "Subject Start", "Subject End", "E-value", "Bit Score"
]
for col in num_cols:
    df[col] = pd.to_numeric(df[col])

df_filtered = df[df["Query ID"] != df["Subject ID"]].copy()

df_filtered["Query Length"] = df_filtered["Query ID"].map(length_lookup)
df_filtered["Subject Length"] = df_filtered["Subject ID"].map(length_lookup)

def compute_vote(row):
    qlen = row["Query Length"]
    return (row["Right Pos"] - row["Left Pos"] + 1) / qlen

df_filtered["Query vote"] = df_filtered.apply(compute_vote, axis=1)

df_comparison = df_filtered[[
    "Query ID", "Subject ID", "Query Length", "Subject Length",
    "Left Pos", "Right Pos", "Subject Start", "Subject End", 
    "Query vote", "Bit Score", "E-value"
]]

df = df_comparison.iloc[:, :-3]

# for CDD input 
cdd_cols = ["Query ID","query length","CDD ID","evalue","q. start","q. end","hit type"]
cdd_rows = []

with open(args.cddi) as fh:
    for ln in fh:
        if ln.startswith('#'):
            continue
        parts = ln.strip().split('\t')

        query_info = parts[1]
        query_id, query_len_str = query_info.rsplit(':', 1)
        query_len = int(query_len_str)

        domain_blocks = parts[2:]

        for block in domain_blocks:
            tokens = block.strip().split('|')
            cdd_id = tokens[0]
            hit_type = tokens[-1]

            if hit_type == "Nohit" and len(tokens) == 2:
                cdd_rows.append([query_id, query_len, cdd_id, '', '', '', 'Nohit'])
                print(f"[INFO] No hit for {query_id} in {cdd_id}")
                continue

            for token in tokens[1:-1]:
                if re.match(r"\d+-\d+:\S+", token):
                    match = re.match(r"(\d+)-(\d+):(\S+)", token)
                    if match:
                        qstart, qend, evalue = match.groups()
                        cdd_rows.append([query_id, query_len, cdd_id, evalue, qstart, qend, hit_type])

cdd_df = pd.DataFrame(cdd_rows, columns=cdd_cols)
cdd_df["evalue"] = pd.to_numeric(cdd_df["evalue"], errors='coerce')
cdd_df[["q. start", "q. end"]] = cdd_df[["q. start", "q. end"]].apply(
    pd.to_numeric, errors='coerce'
)
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

    # Calculate vote
    vote = [set() for _ in range(int(length) + 1)]    

    for row in subset.itertuples():
        subj = row._2
        start = max(1, row._5)
        end = min(length, row._6)
        for i in range(start, end + 1):
            vote[i].add(subj)

    # Identify significant shared regions using vote
    sig_regions = []
    in_region = False
    region_start = None

    for i in range(1, length + 1):
        if (len((vote[i])) / (count - 1))>= args.vote:
            if not in_region:
                region_start = i
                in_region = True
        else:
            if in_region:
                region_end = i - 1
                in_region = False

                # Extend region if it overlaps with a CDD and that CDD has wider boundaries
                for cdd in cdd_subset.itertuples():
                    cdd_start = cdd._5
                    cdd_end = cdd._6

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
            cdd_start = cdd._5
            cdd_end = cdd._6
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
            xmin=row._5, 
            xmax=row._6, 
            colors='green', 
            linewidth=5
        )
        ax.text(length + length*0.1, y_pos, f"{row._3} ({row._7})", va='center', color='black')

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
    ax.set_title(f"Query: {query} (votes >= {args.vote}, n = {count})")

    fig.tight_layout()
    fig.savefig(f"{out_dir}/{query.replace('/', '_')}_sharedRegion.png")
    plt.close()