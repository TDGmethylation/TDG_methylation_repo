import pandas as pd
import numpy as np

def reverse_complement(seq):
    comp_seq = ""
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    for char in seq:
        comp_char = complement[char]
        comp_seq += comp_char

    return comp_seq[::-1]

#Construct the dicticonary
#!!! Should change to hardac, upload input
context_filepath = "/Volumes/Genome_data/TDG_context/TDG_context.txt"

TDG_context = pd.read_table(context_filepath, delimiter="\t",
                            names = ["sequence", "intensity"])
#print(TDG_context)

TDG_dict = {}

#Constant field

#!!! Should change to hardac, upload input
Input_filepath = "/Volumes/Genome_data/TCGA_colon_mutation/Aligned.bed"
Output_basepath = "/Volumes/Genome_data/TCGA_colon_mutation/Result.bed"

Input_format = ["chr", "start", "end", "mut_from", "mut_to", "sequence"]

TCGA_MSI_colon = pd.read_table(Input_filepath, sep = "\t",
                         names = Input_format, index_col = False)

allowed_list = set(["A", "T", "G", "C"])

for index, row in TDG_context.iterrows():
    TDG_dict[row['sequence'][0: 2] + "C" + row['sequence'][3: 7]] = 0

for index, rows in TCGA_MSI_colon.iterrows():

    if((rows['sequence'][6: 7] != "C" and rows['sequence'][6: 7] != "G") or len(rows['sequence']) != 14
            or (allowed_list.issuperset(rows['sequence']) is False)): continue

    if(rows['sequence'][6: 7] == "C"):
        seq = rows['sequence'][4: 11]
        TDG_dict[seq] += 1
        print(seq)

    elif (rows['sequence'][6: 7] == "G"):
        seq = rows['sequence'][2: 9]
        revcomp_seq = reverse_complement(seq)
        TDG_dict[revcomp_seq] += 1

MSI_df = pd.DataFrame.from_dict(TDG_dict, orient = "index")

MSI_df.to_csv(Output_basepath, sep = "\t")