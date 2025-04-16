import pandas as pd
import matplotlib.pyplot as plt
import os

# Load subject match data
df = pd.read_csv("/Users/hazel/Desktop/lab/first_task/data/first_task_output.tsv", sep="\t")
df = df.iloc[:, :-3]

# Load CDD data
cdd_df = pd.read_csv("/Users/hazel/Desktop/lab/first_task/data/cdd_output.tsv", sep="\t")

# Filter for significant CDD hits only
cdd_df = cdd_df[cdd_df["evalue"] < 1e-5]
cdd_df = cdd_df.sort_values(by=["Query ID", "q. start"])

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
df_all.to_csv("/Users/hazel/Desktop/lab/first_task/data/first_task_output_with_reversed.tsv", sep="\t", index=False)

out_dir = "/Users/hazel/Desktop/lab/first_task/query_alignment_graphs"
os.makedirs(out_dir, exist_ok=True)

# for each query
for query in df_all["Query ID"].unique():
    subset = df_all[df_all["Query ID"] == query]
    length = subset["Query Length"].iloc[0]
    
    # Filter CDD annotations for this query
    cdd_subset = cdd_df[cdd_df["Query ID"] == query]

    fig, ax = plt.subplots(figsize=(10, 1.5 + 0.5 * (len(subset) + len(cdd_subset))))

    # plot query sequence at the top
    ax.hlines(
        y=len(subset) + len(cdd_subset) + 1,
        xmin=1,
        xmax=length,
        colors='black',
        linewidth=5,
    )
    # query label
    clean_query_id = query.replace("2.A.123.", "")
    ax.text(length + length*0.1, len(subset) + len(cdd_subset) + 1, clean_query_id, va='center')

    # plot CDD lines on top of subject matches
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
        clean_subject_id = subject.replace("2.A.123.", "")
        ax.text(length + length*0.1, i, clean_subject_id, va='center')

    ax.set_ylim(0, len(subset) + len(cdd_subset) + 2)
    ax.set_xlim(0, length + length*0.35)
    ax.set_xlabel("Query Sequence Position")
    ax.set_yticks([])
    ax.set_title(f"Query: {query}")

    fig.tight_layout()
    fig.savefig(f"{out_dir}/{clean_query_id.replace('/', '_')}_alignment.png")
    plt.close()