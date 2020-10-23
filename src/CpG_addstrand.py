import pandas as pd

#Step1: load alignment data
mut_col = ["chr", "pos_start", "pos_end", "mut_from", "mut_to"]

GT_LAML_KR = pd.read_table("/Volumes/Genome_data/LAML_KR/2nd/CT_GA_filtered.tsv", dtype = object,
                       sep = "\t", names = mut_col, index_col = False)


#Step2: add chr to the chr column (it's just number from the original file) & add strandness

GT_LAML_KR["chr"] = "chr" + GT_LAML_KR["chr"]
strandness = []

for index, row in GT_LAML_KR.iterrows():
    if(row["mut_from"] == "C" and row["mut_to"] == "T"):
        strandness += "+"
    else:
        strandness += "-"

strandnesss = pd.Series(strandness)

GT_LAML_KR = pd.concat([GT_LAML_KR, strandnesss], axis = 1, ignore_index= True)