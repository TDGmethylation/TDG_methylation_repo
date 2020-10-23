import pandas as pd

mut_col = ["chr", "pos_start", "pos_end", "prev_strand",
           "mut_type", "ref", "mut_from", "mut_to", "strand", "pos_end2", "coverage", "methyl", "indicator"]

Merge = pd.read_table("/Volumes/Genome_data/WGBS_K562/data/Merged_rmdup/Merged_rmdup_0_20.tsv",
                           sep="\t", names=mut_col, index_col=False)

col1 = Merge["chr"]
col2 = Merge["pos_start"] - 2
col3 = Merge["pos_start"] + 5
col4 = Merge["mut_from"]
col5 = Merge["mut_to"]
col6 = Merge["strand"]

Input = pd.concat([col1, col2, col3, col4, col5, col6], axis=1, ignore_index=True)

Input.to_csv("/Volumes/Genome_data/WGBS_K562/data/Merged_rmdup/0_20_input.bed", sep="\t", header=False, index=False)