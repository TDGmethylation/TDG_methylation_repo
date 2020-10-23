import numpy as np

import os

ref_list = []

#The global field defines CpG site vs. non-CPG as

CG_high_escore = []
high_escore_total = []
CG_low_escore = []
low_escore_total = []

with open("/Users/gerogehou/Desktop/R/Genomic/GT_parse/GD12_ref/Gli1_v015681_8mers_11111111.txt", "r") as ref:
    head = ref.readline()
    while True:
        current_line = ref.readline()
        if(current_line == ""): break

        current_eachline = current_line.split("\t")
        ref_list.append(current_eachline[0])

marker_index = 0
name_list = []

for filename in os.listdir("/Users/gerogehou/Desktop/R/Genomic/GT_parse/escore"):

    if (str(filename) == ".DS_Store"): continue
    #try:
    #    with open("/Users/gerogehou/Desktop/R/Genomic/GT_parse/escore/" + filename, "r") as escore_file:
    #        for lines in escore_file:
    #            a = 0
    #except UnicodeDecodeError:
    #    print(filename)

    with open("/Users/gerogehou/Desktop/R/Genomic/GT_parse/escore/" + filename, "r", errors = "ignore") as escore_file:

        name_list.append(str(filename))
        CG_high_escore.append(0)
        high_escore_total.append(0)
        CG_low_escore.append(0)
        low_escore_total.append(0)

        ref_index = 0
        #name_list.append(filename)
        #print(filename)

        for lines in escore_file:
            #print(lines, filename)
            #lines = escore_file.readline()
            #print(lines)
            #print(float(lines))
            #lines = lines.decode('utf-8', 'ignore')

            try:
                float(lines)

            except ValueError:
                #print("Not value" + lines)
                continue


            if(float(lines) >= 0.4):
                high_escore_total[marker_index] += 1

                if("CG" in ref_list[ref_index]):
                    CG_high_escore[marker_index] += 1


            else:
                low_escore_total[marker_index] += 1

                if ("CG" in ref_list[ref_index]):
                    CG_low_escore[marker_index] += 1

            ref_index += 1

    marker_index += 1

# print(len(name_list))
# print(name_list)
#
# print(CG_high_escore)
# print(high_escore_total)
# print(CG_low_escore)
# print(low_escore_total)
# print(len(CG_high_escore))
# print(len(high_escore_total))
# print(len(CG_low_escore))
# print(len(low_escore_total))

result_high_CG = []
result_low_CG = []
for i in range(len(CG_low_escore)):
    if high_escore_total[i] != 0 and low_escore_total[i] != 0:
        result_high_CG.append(float(CG_high_escore[i]) / float(high_escore_total[i]))
        result_low_CG.append(float(CG_low_escore[i]) / float(low_escore_total[i]))
    else:
        result_high_CG.append(0.0)
        result_low_CG.append(0.0)


# print(result_high_CG)
# print(result_low_CG)
#
# print(float(sum(CG_high_escore)) / float(sum(high_escore_total)))
# print(float(sum(CG_low_escore)) / float(sum(low_escore_total)))

result_dict = {}
for i in range(len(low_escore_total)):
    print(i)
    result_dict[name_list[i]] = [result_high_CG[i], result_low_CG[i]]

for key in result_dict:
    a = str(result_dict[key][0])
    b = str(result_dict[key][1])
    print(key + "\t" + a + "\t" + b + "\n")


nparray_result = np.empty((667, 3), dtype = object)
r = 0
for key in result_dict:
    nparray_result[r, 0] = key
    nparray_result[r, 1] = result_dict[key][0]
    nparray_result[r, 2] = result_dict[key][1]
    r += 1

np.savetxt("/Users/gerogehou/Desktop/R/Genomic/GT_parse/result/CG.txt", nparray_result, fmt = '%s   %f  %f')