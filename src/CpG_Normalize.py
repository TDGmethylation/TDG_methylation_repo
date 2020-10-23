import pandas as pd

meth_col = ["chr", "pos_start", "pos_end", "strand", "coverage", "methylation_level"]

for i in range(-10, 101, 10):

    InputFilepath = "/Volumes/Genome_data/WGBS_K562/data/bins/Methyl_" + str(i) + "_" + str(i + 10) + ".bed"
    Normalize = pd.read_table(InputFilepath, sep="\t", names=meth_col, index_col=False)

    #print(Normalize)

    chr_col = Normalize["chr"]
    start = Normalize["pos_start"] - 3
    end = Normalize["pos_end"] + 3
    strand = Normalize["strand"]
    void = Normalize["chr"]

    Norm_final = pd.concat([chr_col, start, end, void, void, strand], axis=1, ignore_index=True)

    OutputFilepath = "/Volumes/Genome_data/WGBS_K562/data/alignment_prep/Context_" +\
                     str(i) + "_" + str(i + 10) + ".bed"
    # print(Norm_final)

    Norm_final.to_csv(OutputFilepath, sep="\t", header=False, index=False)

