import numpy as np
import pandas as pd
import math

headers = ["chromosome", "start", "end", "name", "score", "strand", "start_thick",
          "end_thick", "color", "coverage", "methyl_percentage"]

CpG_genome = pd.read_csv("/Users/gerogehou/Desktop/R/genomic/WGBS/intersect_rep.bed", sep = "\t", names = headers)

#print(CpG_genome.head())

methyl_per = CpG_genome['methyl_percentage']

Div_0_10 = pd.DataFrame(columns = headers)
Div_10_20 = pd.DataFrame(columns = headers)
Div_20_30 = pd.DataFrame(columns = headers)
Div_30_40 = pd.DataFrame(columns = headers)
Div_40_50 = pd.DataFrame(columns = headers)
Div_50_60 = pd.DataFrame(columns = headers)
Div_60_70 = pd.DataFrame(columns = headers)
Div_70_80 = pd.DataFrame(columns = headers)
Div_80_90 = pd.DataFrame(columns = headers)
Div_90_100 = pd.DataFrame(columns = headers)


for rows in range(len(methyl_per)):
    if(methyl_per[rows] <= 10):
        Div_0_10 = Div_0_10.append(CpG_genome.iloc[rows])

    elif (10 < methyl_per[rows] <= 20):
        Div_10_20 = Div_10_20.append(CpG_genome.iloc[rows])

    elif (20 < methyl_per[rows] <= 30):
        Div_20_30 = Div_20_30.append(CpG_genome.iloc[rows])

    elif (30 < methyl_per[rows] <= 40):
        Div_30_40 = Div_30_40.append(CpG_genome.iloc[rows])

    elif (40 < methyl_per[rows] <= 50):
        Div_40_50 = Div_40_50.append(CpG_genome.iloc[rows])

    elif (50 < methyl_per[rows] <= 60):
        Div_50_60 = Div_50_60.append(CpG_genome.iloc[rows])

    elif (60 < methyl_per[rows] <= 70):
        Div_60_70 = Div_60_70.append(CpG_genome.iloc[rows])

    elif (70 < methyl_per[rows] <= 80):
        Div_70_80 = Div_70_80.append(CpG_genome.iloc[rows])

    elif (80 < methyl_per[rows] <= 90):
        Div_80_90 = Div_80_90.append(CpG_genome.iloc[rows])

    elif (90 < methyl_per[rows] <= 100):
        Div_90_100 = Div_90_100.append(CpG_genome.iloc[rows])

Div_0_10.to_csv("/Users/gerogehou/Desktop/R/genomic/WGBS/Divided_by_methyl/0_10.csv")
Div_10_20.to_csv("/Users/gerogehou/Desktop/R/genomic/WGBS/Divided_by_methyl/10_20.csv")
Div_20_30.to_csv("/Users/gerogehou/Desktop/R/genomic/WGBS/Divided_by_methyl/20_30.csv")
Div_30_40.to_csv("/Users/gerogehou/Desktop/R/genomic/WGBS/Divided_by_methyl/30_40.csv")
Div_40_50.to_csv("/Users/gerogehou/Desktop/R/genomic/WGBS/Divided_by_methyl/40_50.csv")
Div_50_60.to_csv("/Users/gerogehou/Desktop/R/genomic/WGBS/Divided_by_methyl/50_60.csv")
Div_60_70.to_csv("/Users/gerogehou/Desktop/R/genomic/WGBS/Divided_by_methyl/60_70.csv")
Div_70_80.to_csv("/Users/gerogehou/Desktop/R/genomic/WGBS/Divided_by_methyl/70_80.csv")
Div_80_90.to_csv("/Users/gerogehou/Desktop/R/genomic/WGBS/Divided_by_methyl/80_90.csv")
Div_90_100.to_csv("/Users/gerogehou/Desktop/R/genomic/WGBS/Divided_by_methyl/90_100.csv")
