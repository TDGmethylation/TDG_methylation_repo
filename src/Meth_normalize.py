import pandas as pd
import numpy as np
import math

#Constant field
lower_bound = -10
upper_bound = 0
interval = 10
Back_basepath = "/Volumes/Genome_data/WGBS_K562/data/Background_2nd/Aligned"
Merged_basepath = "/Volumes/Genome_data/WGBS_K562/data/aligned_2nd/Aligned"
Output_basepath = "/Volumes/Genome_data/WGBS_K562/data/Normalization_2nd/Normal"


while lower_bound <= 101:
    Back_filepath = Back_basepath + "_" + str(lower_bound) + "_" +str(upper_bound) + ".bed"
    Merged_filepath = Merged_basepath + "_" + str(lower_bound) + "_" +str(upper_bound) + ".bed"

    Bed_file_inputformat = ["chr", "start", "end", "ignore1", "ignore2", "strand", "sequence"]

    Background = pd.read_table(Back_filepath, dtype = object, sep = "\t",
                             names = Bed_file_inputformat, index_col = False)
    Merged = pd.read_table(Merged_filepath, dtype = object, sep = "\t",
                             names = Bed_file_inputformat, index_col = False)

    #print(Background.shape, Merged.shape)

    Back_positive = Background[Background.strand == "+"].sequence.str.upper()
    Back_negative = Background[Background.strand == "-"].sequence.str.upper()
    Merged_positive = Merged[Merged.strand == "+"].sequence.str.upper()
    Merged_negative = Merged[Merged.strand == "-"].sequence.str.upper()

    #print(Back_positive)

    #Here, background +, we looked at position index = 3
    #background -, we looked at position index = 4
    #merged +, we looked at position index = 3
    #merged -, we looked at position index = 4

    Back_nucleotide_count = {"A":0, "C":0, "G":0, "T":0}
    Merged_nucleotide_count = {"A": 0, "C": 0, "G": 0, "T": 0}

    #print(type(Back_positive))

    for sequence in Back_positive:
        Back_nucleotide_count[sequence[3]] += 1

    for sequence in Back_negative:
        Back_nucleotide_count[sequence[4]] += 1

    for sequence in Merged_positive:
        Merged_nucleotide_count[sequence[3]] += 1

    for sequence in Merged_negative:
        Merged_nucleotide_count[sequence[4]] += 1

    Back_count = pd.Series(Back_nucleotide_count, name = "Background", dtype = float)
    Merged_count = pd.Series(Merged_nucleotide_count, name = "Merged", dtype = float)
    Divided = Merged_count/Back_count
    Divided.rename("Ratio")
    print(Divided)
    #print(Back_count)
    #print(Merged_count)

    Result = pd.concat([Back_count, Merged_count, Divided], axis = 1)
    Result_append = Result
    #print(Result)

    Output = Output_basepath + "_" + str(lower_bound) + "_" +str(upper_bound) + ".tsv"

    Result.to_csv(Output, sep="\t", header=True, index=True)

    #print(Back_nucleotide_count, Merged_nucleotide_count, lower_bound, upper_bound)

    lower_bound += interval
    upper_bound += interval




