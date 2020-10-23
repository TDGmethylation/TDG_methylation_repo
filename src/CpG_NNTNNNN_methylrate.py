import pandas as pd
import numpy as np

#Construct the dicticonary
#!!! Should change to hardac, upload input
context_filepath = "/Volumes/Genome_data/TDG_context/TDG_context.txt"

TDG_context = pd.read_table(context_filepath, delimiter="\t",
                            names = ["sequence", "intensity"])
#print(TDG_context)

TDG_dict_Back = {}
TDG_dict_MetRatio = {}

#Constant field

#!!! Should change to hardac, upload input
Input_filepath = "/Volumes/Genome_data/WGBS_NeuronStemCell/Aligned_NSC_cov10.bed"
Output_basepath = "/Volumes/Genome_data/WGBS_NeuronStemCell/Result/Result_"

Back_all = pd.DataFrame()
MetRatio_all = pd.DataFrame()

Input_format = ["chr", "start", "end", "coverage", "methyl_ratio", "strand", "sequence"]


hg19_aligned = pd.read_table(Input_filepath, sep = "\t",
                         names = Input_format, index_col = False)

# print(hg19_aligned.shape)

hg19_positive = hg19_aligned[hg19_aligned.strand == "+"]
hg19_positive['sequence'] = hg19_positive['sequence'].str.upper()

hg19_negative = hg19_aligned[hg19_aligned.strand == "-"]
hg19_negative['sequence'] = hg19_negative['sequence'].str.upper()

# hg19_positive = hg19_aligned[hg19_aligned.strand == "+"].sequence.str.upper()
# hg19_negative = hg19_aligned[hg19_aligned.strand == "-"].sequence.str.upper()

# print(hg19_positive.shape, hg19_negative.shape)
# print(hg19_positive)

allowed_list = set(["A", "T", "G", "C"])


#for coverage_cutoff in range(0, 50, 5):
for coverage_cutoff in range(10, 50, 5):

    for index, row in TDG_context.iterrows():
        TDG_dict_Back[row['sequence'][0: 2] + "C" + row['sequence'][3: 7]] = 0
        TDG_dict_MetRatio[row['sequence'][0: 2] + "C" + row['sequence'][3: 7]] = 0

    for index, rows in hg19_positive.iterrows():

        #Filter reads that are not methylated
        if (float(rows['methyl_ratio']) <= 0.1): continue

        if(rows["coverage"] < coverage_cutoff): continue
        #print(rows["coverage"], coverage_cutoff)
        if(rows['sequence'][6: 7] != "C" or len(rows['sequence']) != 14 or
                allowed_list.issuperset(rows['sequence']) is False): continue

        #print(rows['sequence'][4: 11])
        TDG_dict_Back[rows['sequence'][4: 11]] += 1
        #print(rows['methyl_ratio'])

        #Previously, I put percentage 20% as 20, so there is a float(rows['methyl_ratio'])/100.0
        #there. But this this H1 cell line, I put 0.2 as 20. so I remove the 100.0.
        #Just to be aware.

        TDG_dict_MetRatio[rows['sequence'][4: 11]] += float(rows['methyl_ratio'])

        # Indicator of progress
        #if(index % 1000 == 0): print("Reach %d at positive" % index)

    for index, rows in hg19_negative.iterrows():

        #same as above
        if (float(rows['methyl_ratio']) <= 0.1): continue


        if (rows["coverage"] < coverage_cutoff): continue
        if(rows['sequence'][7: 8] != "C" or len(rows['sequence']) != 14 or
                allowed_list.issuperset(rows['sequence']) is False): continue

        #print(sequence[2: 9])
        TDG_dict_Back[rows['sequence'][5: 12]] += 1
        TDG_dict_MetRatio[rows['sequence'][5: 12]] += float(rows['methyl_ratio'])


        #Indicator of progress
        #if(index % 1000 == 0): print("Reach %d at negative" % index)

    #print(TDG_dict_Back, TDG_dict_Merge)

    # Merged_Output = Output_basepath + "Merged_" + str(lower_bound) + "_" + str(upper_bound) + ".bed"
    # Back_Output = Output_basepath + "Background_" + str(lower_bound) + "_" + str(upper_bound) + ".bed"

    MetRatio_Output_df = pd.DataFrame.from_dict(TDG_dict_MetRatio, orient = "index")

    Back_Output_df = pd.DataFrame.from_dict(TDG_dict_Back, orient = "index")

    #print(Merged_Output_df)

    MetRatio_Output = Output_basepath + "coverage_" + str(coverage_cutoff) + "_MetRatio.bed"
    Back_Output = Output_basepath + "coverage_" + str(coverage_cutoff) + "_Background.bed"

    MetRatio_Output_df.to_csv(MetRatio_Output, sep = "\t")
    Back_Output_df.to_csv(Back_Output, sep = "\t")