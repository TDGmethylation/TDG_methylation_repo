import numpy as np
import io
import toolkit
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import random

#loop_generator generates a list of all n-mer loops recursivly. All 5'-3'.

design_context = toolkit.Sequence_dict_initializer_list(6, "T", 2)

Loop_context_v1 = []
Loop_context_v1_name = []

loop_context = []

hybrid_clamp = "CGC"
design_precontext = "AGCTC"
design_postcontext = "CT"
loop_clamp = "CGC"
loop = "GGAA"
rev_loop_clamp = "GCG"
rev_design_postcontext = "AG"
rev_design_precontext = "GAGCT"
rev_hybrid_clamp = "GCG"
ssDNA = "GTAGACATTAAAGATT"

# rev_test = Seq(ssDNA, IUPAC.unambiguous_dna)
# reverseda = rev_test.complement()
# print(reverseda)

count = 0

#part 1: design the NNXNNNN main part. N = 40960

for sequence in design_context:
    Seqformat_sequence = Seq(sequence, IUPAC.unambiguous_dna)
    rev_sequence_wc = Seqformat_sequence.reverse_complement()
    rev_sequence_mm = rev_sequence_wc[: 4] + "G" + rev_sequence_wc[-2 :]
    #print(rev_sequence_mm, sequence)

    #the sixty mer complete sequence from 5'-3'. sequence varibale is the desisn NNXNNNN part.
    sixty_mer = hybrid_clamp + design_precontext + sequence + design_postcontext + loop_clamp \
    + loop + rev_loop_clamp + rev_design_postcontext + rev_sequence_mm + rev_design_precontext \
    + rev_hybrid_clamp + ssDNA
    for rep in range(10):
        loop_context.append(sixty_mer)

        #Name = LoopContext1_NNXNNNN_r1-10
        Loop_context_v1_name.append("LoopContext1_" + str(sequence) + "_r" + str(rep + 1))

        count += 1
        #print(sixty_mer)

#print(count)

#part2.1 WT sequence control 173 x 10 = 18200

sequence_WTcontrol = []

for i in range(182):
    ran_index = random.randint(0, count / 10)
    rand_seq_mm = loop_context[ran_index * 10]
    rand_seq_wt = rand_seq_mm[:10] + "C" + rand_seq_mm[11:]
    for rep in range(10):
        #print(rand_seq_wt)
        sequence_WTcontrol.append(rand_seq_wt)

        # Name = LoopContext1_wt_i_r1-10
        Loop_context_v1_name.append("LoopContext1_wt" + str(i + 1) + "_r" + str(rep + 1))

#part2.2 MutS binding ssDNA control

sequence_MutS = []

hybrid_clamp = "CGC"
x_part = ""
MutS_mm = "TGTCTGACT"
MutS_wt = "TGCCTGACT"
loop_clamp = "CGC"
loop = "GGAA"
rev_loop_clamp = "GCG"
rev_MutS = "AGTCAGGCA"
rev_x_part = ""
rev_hybrid_clamp = "GCG"
ssDNA = ""

#this ssDNA is from Harshit, but at 5' I add TT, no influence by mfold
ssDNA_RefList = "TTAACGATGAGTAGACATTAAAGATT"

# x_part_RefList = "AGCTCAAAAAAAA"
# rev_x_part_RefList = "TTTTTTTTGAGCT"

x_part_RefList = "AGCTCATATATAT"
rev_x_part_RefList = "ATATATATGAGCT"

#x is the length of x part ranging from 1 to 13.
for x in range(14):
    x_part = x_part_RefList[: x]
    rev_x_part = rev_x_part_RefList[13 - x :]
    #print(x_part + "\n" + rev_x_part)
    ssDNA = ssDNA_RefList[2 * x :]
    #print(ssDNA)
    sixty_mer_mm = hybrid_clamp + x_part + MutS_mm + loop_clamp + loop + rev_loop_clamp + \
                rev_MutS + rev_x_part + rev_hybrid_clamp + ssDNA

    sixty_mer_wt = hybrid_clamp + x_part + MutS_wt + loop_clamp + loop + rev_loop_clamp + \
                rev_MutS + rev_x_part + rev_hybrid_clamp + ssDNA

    #print(sixty_mer_mm + "\n" + sixty_mer_wt)

    for rep in range(10):
        sequence_MutS.append(sixty_mer_mm)

        # Name = LoopContext1_MutSmm_ss_#_r1-10
        Loop_context_v1_name.append("LoopContext1_MutSmm_ss_" + str(len(ssDNA)) + "_r" + str(rep + 1))

    for rep in range(10):
        sequence_MutS.append(sixty_mer_wt)

        # Name = LoopContext1_MutSwt_ss_#_r1-10
        Loop_context_v1_name.append("LoopContext1_MutSwt_ss_" + str(len(ssDNA)) + "_r" + str(rep + 1))



#part2.3 TDG & hybridization: different lengths of ssDNA and HE.

sequence_HE = []

hybrid_clamp = "CGC"
y_part = ""
TDG_mm = "GGTGCCTCT"
TDG_wt = "GGCGCCTCT"
loop_clamp = "CGC"
loop = "GGAA"
rev_loop_clamp = "GCG"
rev_TDG = "AGAGGCGCC"
rev_y_part = ""
rev_hybrid_clamp = "GCG"
ssDNA_HE = ""

# this ssDNA is from Harshit, but at 5' I add TT, no influence by mfold
ssDNA_RefList_HE = "TTAACGATGAGTAGACATTAAAGATT"
#
# y_part_RefList_HE = "AGCTCAAAAAAAA"
# rev_y_part_RefList_HE = "TTTTTTTTGAGCT"

y_part_RefList_HE = "AGCTCATATATAT"
rev_y_part_RefList_HE = "ATATATATGAGCT"

# x is the length of x part ranging from 1 to 13.
for y in range(14):
    y_part = y_part_RefList_HE[: y]
    rev_y_part = rev_y_part_RefList_HE[13 - y:]
    # print(y_part + "\n" + rev_y_part)
    ssDNA_HE = ssDNA_RefList_HE[2 * y:]
    # print(ssDNA)
    sixty_mer_mm_HE = hybrid_clamp + y_part + TDG_mm + loop_clamp + loop + rev_loop_clamp + \
                   rev_TDG + rev_y_part + rev_hybrid_clamp + ssDNA_HE

    sixty_mer_wt_HE = hybrid_clamp + y_part + TDG_wt + loop_clamp + loop + rev_loop_clamp + \
                   rev_TDG + rev_y_part + rev_hybrid_clamp + ssDNA_HE

    #print(sixty_mer_mm_HE + "\n" + sixty_mer_wt_HE)

    for rep in range(10):
        sequence_HE.append(sixty_mer_mm_HE)
        # Name = LoopContext1_HEmm_ss#_r1-10
        Loop_context_v1_name.append("LoopContext1_HEmm_ss_" + str(len(ssDNA_HE)) + "_r" + str(rep + 1))

    for rep in range(10):
        sequence_HE.append(sixty_mer_wt_HE)

        # Name = LoopContext1_HEwt_ss#_r1-10
        Loop_context_v1_name.append("LoopContext1_HEwt_ss_" + str(len(ssDNA_HE)) + "_r" + str(rep + 1))

#2.4 controls for loops

loop_control = []

preloop = "CGCAGCTCGGTGCCTCTCGC"
preloop_wt = "CGCAGCTCGGCGCCTCTCGC"
loopcontrol_loop = ["GGAA", "CGTA", "TCCT", "CTTT", "TGTC", "GGTC", "CGTG", "ATGG"]
postloop = "GCGAGAGGCGCCGAGCTGCG"
ssDNA_loop = "GTAGACATTAAAGATT"

for loop_comp in loopcontrol_loop:
    for i in range(10):
        loop_sixtymer = preloop + loop_comp + postloop + ssDNA_loop
        loop_control.append(loop_sixtymer)

        # Name = Loop_context1_loop_XXXX_r1-10
        Loop_context_v1_name.append("Loop_context1_loop_" + loop_comp + "_mm_r" + str(i + 1))

        #print(loop_sixtymer)

    for i in range(10):
        loop_sixtymer_wt = preloop_wt + loop_comp + postloop + ssDNA_loop
        loop_control.append(loop_sixtymer_wt)

        # Name = Loop_context1_loop_XXXX_r1-10
        Loop_context_v1_name.append("Loop_context1_loop_" + loop_comp + "_wt_r" + str(i + 1))

        # print(loop_sixtymer)

#2.5 T-G position control for TDG !!!NEED TO ELIMINATE SEVERAL PROBES
TG_pos = []
TDG_binding = "TGTGC"
rev_TDG_binding = "GCGCA"
TDG_binding_wt = "TGCGC"

hybrid_clamp = "CGC"
insert_seq = "GCACGTCCT"
loop_clamp = "CGC"
loop = "GGAA"
rev_loop_clamp = "GCG"
rev_insert_seq = "AGGACGTGC"
rev_hybrid_clamp = "GCG"

ssDNA_pos = "GTAGACATTAAAGATT"

for index in range(10):
    inserted_seq = insert_seq[: index] + TDG_binding + insert_seq[index :]
    rev_inserted_seq = rev_insert_seq[index :] + rev_TDG_binding + rev_insert_seq[: index]
    inserted_seq_wt = insert_seq[: index] + TDG_binding_wt + insert_seq[index:]
    #print(inserted_seq + '\n' + rev_inserted_seq)

    sixtymer_pos_mm = hybrid_clamp + inserted_seq + loop_clamp + loop + rev_loop_clamp + rev_inserted_seq \
    + rev_hybrid_clamp + ssDNA_pos
    sixtymer_pos_wt = hybrid_clamp + inserted_seq_wt + loop_clamp + loop + rev_loop_clamp + rev_inserted_seq \
    + rev_hybrid_clamp + ssDNA_pos

    for rep in range(10):
        TG_pos.append(sixtymer_pos_mm)

        # Name = LoopContext1_TGpos_"# of nt to the loop"_r1-10
        Loop_context_v1_name.append("LoopContext1_TGPosmm_" + str(12 - index)  + "_r" + str(rep + 1))

    for rep in range(10):
        TG_pos.append(sixtymer_pos_wt)

        # Name = LoopContext1_TGpos_"# of nt to the loop"_r1-10
        Loop_context_v1_name.append("LoopContext1_TGPoswt_" + str(12 - index)  + "_r" + str(rep + 1))


#2.6 shrinking ssDNA probes and complete ssDNA probes

sequence_MutS_shrinkssDNA = []

Constant_MutS_mm = "CGCAGCTCTGTCTGACTCGCGGAAGCGAGTCAGGCAGAGCTGCG"
Constant_MutS_wt = "CGCAGCTCTGCCTGACTCGCGGAAGCGAGTCAGGCAGAGCTGCG"
Constant_ssDNA = "GTAGACATTAAAGATT"

for i in range(16, -1, -4):
    partial_ssDNA = Constant_ssDNA[i: 16]
    #print(partial_ssDNA + "end")
    Constant_probes_mm = Constant_MutS_mm + partial_ssDNA
    Constant_probes_wt = Constant_MutS_wt + partial_ssDNA
    ssDNA_length = 16 - i

    for rep in range(10):
        sequence_MutS_shrinkssDNA.append(Constant_probes_mm)
        Loop_context_v1_name.append("LoopContext1_shrink_ssDNA_" + str(ssDNA_length) + "_mm_r" + str(rep + 1))

    for rep in range(10):
        sequence_MutS_shrinkssDNA.append(Constant_probes_wt)
        Loop_context_v1_name.append("LoopContext1_shrink_ssDNA_" + str(ssDNA_length) + "_wt_r" + str(rep + 1))

# for i in range(16, -1, -4):
#     partial_ssDNA = Constant_ssDNA[i: 16]
#     ssDNA_length = 16 - i
#     for rep in range(10):
#         sequence_MutS_shrinkssDNA.append(partial_ssDNA)
#         Loop_context_v1_name.append("LoopContext1_shrink_ssDNA_" + str(ssDNA_length) + "_only_r" + str(rep + 1))

# #2.7 blank
# blank = ["", "", ""]
#
# # Name = LoopContext1_blank_r1-10
# Loop_context_v1_name.append("LoopContext1_blank_r1")
# Loop_context_v1_name.append("LoopContext1_blank_r2")
# Loop_context_v1_name.append("LoopContext1_blank_r3")


#Summary
Loop_context_v1.extend(loop_context)
Loop_context_v1.extend(sequence_WTcontrol)
Loop_context_v1.extend(sequence_MutS)
Loop_context_v1.extend(sequence_HE)
Loop_context_v1.extend(loop_control)
Loop_context_v1.extend(TG_pos)
Loop_context_v1.extend(sequence_MutS_shrinkssDNA)
# Loop_context_v1.extend(blank)

print(len(Loop_context_v1))


with open("/Users/gerogehou/Desktop/R/Array_Design/Loop_context_v1/LoopContext_Seq.txt", "w") as Seq_file:
    for items in Loop_context_v1:

        seq_result = Seq._get_seq_str_and_check_alphabet(items, items)
        Seq_file.write(seq_result + "\n")

with open("/Users/gerogehou/Desktop/R/Array_Design/Loop_context_v1/LoopContext_Name.txt", "w") as Name_file:
    for items in Loop_context_v1_name:
        #print(items)
        name_result = Seq._get_seq_str_and_check_alphabet(items, items)
        Name_file.write(name_result + "\n")
