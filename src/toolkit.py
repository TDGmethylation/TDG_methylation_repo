

def Sequence_dict_initializer_map(length, fixed, position):
    def N_generator(randomized_length):
        nucleotide_tuple = ("A", "T", "G", "C")
        # Base case
        if randomized_length == 1:
            return ["A", "T", "G", "C"]
        # Recursive call
        result1 = N_generator(randomized_length - 1)
        temp_result = []
        for seq_index in range(len(result1)):
            for i in range(4):
                temp_result.append(nucleotide_tuple[i] + result1[seq_index])
        return temp_result
    two_mer_list = N_generator(length)
    six_mer_map_append = {}

    mark = 0
    for twomer in two_mer_list:
        twomer = twomer[:position] + fixed + twomer[position:]
        #print(twomer)
        six_mer_map_append[twomer] = 0
    return six_mer_map_append

def Sequence_dict_initializer_list(length, fixed, position):
    def N_generator(randomized_length):
        nucleotide_tuple = ("A", "T", "G", "C")
        # Base case
        if randomized_length == 1:
            return ["A", "T", "G", "C"]
        # Recursive call
        result1 = N_generator(randomized_length - 1)
        temp_result = []
        for seq_index in range(len(result1)):
            for i in range(4):
                temp_result.append(nucleotide_tuple[i] + result1[seq_index])
        return temp_result
    two_mer_list = N_generator(length)
    six_mer_list_append = []
    mark = 0
    for twomer in two_mer_list:
        twomer = twomer[:position] + fixed + twomer[position:]
        #print(twomer)
        six_mer_list_append.append(twomer)
    return six_mer_list_append

#print(Sequence_dict_initializer_map(2, "T", 0))