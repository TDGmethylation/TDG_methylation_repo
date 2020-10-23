import json
import os
import pandas as pd

root = "/Volumes/Genome_data/TCGA_colon_mutation/"
summary_df = pd.DataFrame()

for filenames in os.listdir(root):
    if(filenames.endswith(".json")):
        with open(root + filenames) as f:
            TCGA = json.load(f)
            for mutation_type in TCGA:
                notation = mutation_type["genomic_dna_change"]
                if(">" in notation):
                    splited_notation = notation.split(":")
                    chr = splited_notation[0]
                    splited_second = splited_notation[1].split(">")
                    index = int(splited_second[0][2: -1])
                    mut_from = splited_second[0][-1: ]
                    mut_to = splited_second[1]
                    pre_series = [chr, index, index+1, mut_from, mut_to]
                    #print(pre_series)
                    bed_series = pd.Series(pre_series)
                    summary_df = pd.concat([summary_df, bed_series], axis=1)
                    #print(chr)

summary_df = summary_df.T
#summary_df.col = ["chr", "start", "end", "mut_from", "mut_to"]
path = "/Volumes/Genome_data/TCGA_colon_mutation/Aggregrated_mutation.bed"
summary_df.to_csv(path, sep='\t', header=False, index=False)
