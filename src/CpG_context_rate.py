import pandas as pd
import numpy as np

#Construct the dicticonary
context_filepath = "/Volumes/Genome_data/TDG_context/TDG_context.txt"

TDG_context = pd.read_table(context_filepath, dtype = object, delimiter="\t",
                            names = ["sequence", "intensity"])
#print(TDG_context)

TDG_dict_Back = {}
TDG_dict_Merge = {}

#Constant field
lower_bound = 0
upper_bound = 10
interval = 10
Back_basepath = "/Volumes/Genome_data/WGBS_K562/raw/Unioned/CN/CN_category/Cat"
Merged_basepath = "/Volumes/Genome_data/WGBS_K562/data_CpN/2_Categorize/Cat"
Output_basepath = "/Volumes/Genome_data/WGBS_K562/data_CpN/5_mut_ratio/"

Back_all = pd.DataFrame()
Merged_all = pd.DataFrame()


while lower_bound < 100:
    print("Current in step" + str(lower_bound))
    for index, row in TDG_context.iterrows():
        TDG_dict_Back[row['sequence'][0: 2] + "C" + row['sequence'][3: 7]] = 0
        TDG_dict_Merge[row['sequence'][0: 2] + "C" + row['sequence'][3: 7]] = 0

    Back_filepath = Back_basepath + "_" + str(lower_bound) + "_" +str(upper_bound) + ".bed"
    Merged_filepath = Merged_basepath + "_" + str(lower_bound) + "_" +str(upper_bound) + ".tsv"

    Back_file_inputformat = ["chr", "start", "end", "ignore1", "ignore2", "strand", "sequence"]

    Merged_file_inputformat = ["chr", "start", "end", "ig1", "ig2", "strand", "index",
                               "ig3", "ig4", "ig5", "sequence"]

    Background = pd.read_table(Back_filepath, dtype = object, sep = "\t",
                             names = Back_file_inputformat, index_col = False)
    Merged = pd.read_table(Merged_filepath, dtype = object, sep = "\t",
                             names = Merged_file_inputformat, index_col = False)

    print(Background.shape, Merged.shape)

    Back_positive = Background[Background.strand == "+"].sequence.str.upper()
    Back_negative = Background[Background.strand == "-"].sequence.str.upper()
    Merged_positive = Merged[Merged.strand == "+"].sequence.str.upper()
    Merged_negative = Merged[Merged.strand == "-"].sequence.str.upper()

    #print(Back_positive)

    allowed_list = {"A", "C", "G", "T"}

    for sequence in Back_positive:
        if(sequence[6: 7] != "C" or len(sequence) != 14 or
                allowed_list.issuperset(sequence) is False): continue

        #print(sequence[2: 9])
        TDG_dict_Back[sequence[4: 11]] += 1
        #print(TDG_dict_Back[sequence[2: 9]])

    for sequence in Back_negative:
        if (sequence[7: 8] != "C" or len(sequence) != 14 or
                allowed_list.issuperset(sequence) is False): continue
        TDG_dict_Back[sequence[5: 12]] += 1

    for sequence in Merged_positive:
        if (sequence[6: 7] != "C" or len(sequence) != 14 or
                allowed_list.issuperset(sequence) is False): continue
        TDG_dict_Merge[sequence[4: 11]] += 1

    for sequence in Merged_negative:
        if (sequence[7: 8] != "C" or len(sequence) != 14 or
                allowed_list.issuperset(sequence) is False): continue
        TDG_dict_Merge[sequence[5: 12]] += 1

    #print(TDG_dict_Back, TDG_dict_Merge)

    # Merged_Output = Output_basepath + "Merged_" + str(lower_bound) + "_" + str(upper_bound) + ".bed"
    # Back_Output = Output_basepath + "Background_" + str(lower_bound) + "_" + str(upper_bound) + ".bed"

    Segments = str(lower_bound) + "_" +str(upper_bound)

    Merged_Output_df = pd.DataFrame.from_dict(TDG_dict_Merge, orient = "index")
    Merged_Output_df.columns = [Segments]
    #print(Merged_Output_df)
    Merged_all = pd.concat([Merged_all, Merged_Output_df], axis = 1)


    Back_Output_df = pd.DataFrame.from_dict(TDG_dict_Back, orient = "index")
    Back_Output_df.columns = [Segments]
    Back_all = pd.concat([Back_all, Back_Output_df], axis=1)

    #print(Merged_Output_df)


    lower_bound += interval
    upper_bound += interval

Merged_All_Output = Output_basepath + "Intersect_All.bed"
Back_All_Output = Output_basepath + "Background_All.bed"
Merged_all.to_csv(Merged_All_Output, sep = "\t")
Back_all.to_csv(Back_All_Output, sep="\t")