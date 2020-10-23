import pandas as pd
import numpy as np

Intersected_filepath = "/Volumes/Genome_data/WGBS_K562/data_CpN_byPosition/All_positional.csv"
Outpath = "/Volumes/Genome_data/WGBS_K562/data_CpN_byPosition/Intersected_removeNaN_devided10.tsv"

Intersected = pd.read_csv(Intersected_filepath, sep=",")

print(Intersected)

Intersect_namelist = ["0_20", "20_40", "40_60", "60_80", "80_100"]

Output_intersected = pd.DataFrame(columns=Intersect_namelist)

for index, row in Intersected.iterrows():
    for meth in range(0, 99, 20):
        header = str(meth) + "_" + str(meth + 20)
        #print(header)
        appendable_series = pd.DataFrame(row.Intensity * np.ones(int(row[header]/10)), columns=[header])
        #print(appendable_series)
        Output_intersected = pd.concat([Output_intersected, appendable_series], axis=0,
                                       ignore_index=True)
        #print(Output_intersected)
    print(index)

remove_NaN = (Output_intersected.apply(lambda x: pd.Series(x.dropna().values)))

remove_NaN.to_csv(Outpath, sep="\t", header=False, index=False)
