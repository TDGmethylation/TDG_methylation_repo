import pandas as pd


for lo in range(-10, 101, 10):
    input_lower = lo
    input_upper = lo + 10

    # Step1: load alignment data
    mut_col = ["chr", "pos_start", "pos_end", "prev_strand",
               "mut_type", "ref", "mut_from", "mut_to", "strand"]

    GT_LAML_KR = pd.read_table("/Volumes/Genome_data/LAML_KR/GT_CML_plusnote_sorted.tsv", dtype=object,
                               sep="\t", names=mut_col, index_col=False)

    #print(GT_LAML_KR.head(10))

    methyl_col = ["chr", "pos_start_minus_1", "pos_start", "strand", "coverage", "meth_level_falsehere"]
    methyl = pd.read_table("/Volumes/Genome_data/WGBS_K562/data/liftover_20bp/Shifted_" + str(input_lower) + "_" +
                           str(input_upper) + ".bed", dtype=object, sep="\t", names=methyl_col, index_col=False)
    # print(methyl.head(10))

    # Step3: Merge two files with chr and pos_start

    Merged = pd.merge(GT_LAML_KR, methyl, on=["chr", "pos_start", "strand"], how='inner', indicator=True)
    print(Merged.head(10))

    # Merged.to_csv("/Volumes/Genome_data/WGBS/" + str(input_lower)
    #              + "_" + str(input_upper) + "_merged.csv")
    # chr_series = Merged["chr"]
    # start = Merged["pos_start"]
    # start = start - 2
    # end = Merged["pos_start"]
    # end = end + 5
    # strand = Merged["strand"]
    # bed_input = pd.concat([chr_series, start, end], axis=1, ignore_index=True)
    # #print(bed_input.head(10))

    Merged.to_csv("/Volumes/Genome_data/WGBS_K562/data/Merge_2nd/MergeUpdate" + "_" + str(input_lower)
                  + "_" + str(input_upper) + ".tsv", sep="\t", header=False, index=False)
