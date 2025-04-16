import pandas as pd
import re

file_path = "/Users/hazel/Desktop/lab/first_task/data/ssearch.out" 

query_lengths = {} 
filtered_lines = [] 

with open(file_path, "r") as file:
    for line in file:
        # store query length
        match = re.match(r"# Query: (\S+)  - (\d+) aa", line)
        if match:
            query_id = match.group(1)
            query_lengths[query_id] = int(match.group(2))

        # store only non-header lines
        elif not line.startswith("#"):
            filtered_lines.append(line)

columns = ["Query ID", "Subject ID", "% Identity", "Alignment Length", "Mismatches", "Gap Opens", 
           "Left Pos", "Right Pos", "Subject Start", "Subject End", "E-value", "Bit Score"]

data = [line.split("\t") for line in filtered_lines]
df = pd.DataFrame(data, columns=columns)

# Convert to num
num_cols = ["% Identity", "Alignment Length", "Mismatches", "Gap Opens", 
            "Left Pos", "Right Pos", "Subject Start", "Subject End", "E-value", "Bit Score"]
for col in num_cols:
    df[col] = pd.to_numeric(df[col]) 

# Remove self comparisons
df_filtered = df[df["Query ID"] != df["Subject ID"]].copy()

# Compute Qcov
df_filtered["Query Coverage"] = df_filtered.apply(
    lambda row: (row["Right Pos"] - row["Left Pos"] + 1) / query_lengths.get(row["Query ID"]), axis=1
)

df_filtered["Query Length"] = df_filtered["Query ID"].map(query_lengths)
df_filtered["Subject Length"] = df_filtered["Subject ID"].map(query_lengths)  # ← this is new

df_comparison = df_filtered[["Query ID", "Subject ID", "Query Length", "Subject Length",
                             "Left Pos", "Right Pos", "Subject Start", "Subject End", 
                             "Query Coverage", "Bit Score", "E-value"]]

df_comparison.to_csv("first_task_output.tsv", sep="\t", index = False)
