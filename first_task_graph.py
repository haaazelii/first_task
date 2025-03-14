import pandas as pd
import matplotlib.pyplot as plt
import os

df = pd.read_csv("first_task_output.tsv", sep="\t")
out_dir = "query_alignment_graphs"
os.makedirs(out_dir, exist_ok=True)

# for each query
for query in df["Query ID"].unique():
    # subset for current query
    subset = df[df["Query ID"] == query]
    length = subset["Query Length"].iloc[0]

    fig, ax = plt.subplots(figsize=(10, 1.5 + 0.5 * len(subset)))

    # plot query sequence at the top
    ax.hlines(
        y=len(subset) + 1,
        xmin=1,
        xmax=length,
        colors='black',
        linewidth=5,
    )
    # query label
    clean_query_id = query.replace("2.A.123.", "")
    ax.text(length + length*0.1, len(subset) + 1, clean_query_id, va='center')

    # plot subjects
    for i, row in enumerate(reversed(list(subset.itertuples())), 1):
        left_pos = row._4 
        right_pos = row._5  
        subject = row._2 
        ax.hlines(
            y=i, 
            xmin=left_pos, 
            xmax=right_pos, 
            colors='blue', 
            linewidth=5
        )
        # subject label
        clean_subject_id = subject.replace("2.A.123.", "")
        ax.text(length + length*0.1, i, clean_subject_id, va='center')

    ax.set_ylim(0, len(subset) + 2)
    ax.set_xlim(0, length + length*0.35)
    ax.set_xlabel("Query Sequence Position")
    ax.set_yticks([])
    ax.set_title(f"Query: {query}")

    fig.tight_layout()
    fig.savefig(f"{out_dir}/{clean_query_id.replace('/', '_')}_alignment.png")
    plt.close()
