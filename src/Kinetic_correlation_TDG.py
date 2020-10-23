import pandas as pd
import math

filepath = "/Users/gerogehou/Desktop/R/Correlation_with_papers/correlation_input5.csv"
list = ["sequence", "kmax", "intensity"]
TDG_table = pd.read_csv(filepath, sep=',', names=list)

output_t0 = []
ratio = {}

#print(TDG_table)

#Here you can change the loop param to select columns. a in range (0, 8), 前8行。

for a in range(0, 5):
    for b in range(a + 1, 6):
        p = TDG_table.at[a, "intensity"] / TDG_table.at[b, "intensity"]
        ka = TDG_table.at[a, "kmax"]
        kb = TDG_table.at[b , "kmax"]
        # t0 = 2 * (ka - kb * p) / (ka * ka - p * kb * kb)
        # output_t0.append(t0)

        t0 = 0.05

        exp = (1 - math.pow(math.e, -1 * ka * t0)) / (1 - math.pow(math.e, -1 * kb * t0))

        ratio[p] = exp

# print(output_t0)
print(ratio)