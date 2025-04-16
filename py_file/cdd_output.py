import csv

input_file = "/Users/hazel/Desktop/lab/first_task/data/cdd.out"
output_file = "/Users/hazel/Desktop/lab/first_task/data/cdd_output.tsv"

columns = [
    "Query ID", "query length", "CDD ID", "CDD length", "evalue", "bit score",
    "q. start", "q. end", "% query coverage per hsp", "s. start", "s. end"
]

with open(input_file) as infile, open(output_file, "w", newline="") as outfile:
    writer = csv.writer(outfile, delimiter="\t")
    writer.writerow(columns)

    for line in infile:
        if line.startswith("#"):
            continue

        fields = line.strip().split("\t")
        if len(fields) >= 11:
            writer.writerow(fields[:11])