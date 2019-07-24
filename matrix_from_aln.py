import numpy


def create_matrix_from_aln(filename):
    with open(filename) as file_handle:
        lines = {}
        for index, line in enumerate(file_handle):
            tokenized_line = line.split()
            if index > 0 and len(tokenized_line) > 1:
                if tokenized_line[0] not in lines:
                    lines[tokenized_line[0]] = list(tokenized_line[1].lower())
                else:
                    lines[tokenized_line[0]] += list(tokenized_line[1].lower())
        in_array = []
        for key, line in lines.items():
            print(line)
            in_array.append(line)

        return numpy.array(in_array, numpy.unicode_)


def load_sequences_from_fasta(filename):
    lines = []
    with open(filename) as file_handle:
        for line in file_handle:
            if line[0] != '>':
                lines.append(line[0:-1])

    return lines


def matrix_from_fasta(filename):
    lines = load_sequences_from_fasta(filename)
    list_of_lists = []
    for line in lines:
        converted = list(line)
        list_of_lists.append(converted)
        print(converted)
    return list_of_lists


def matrix_from_exa(filename):
    list_of_list = []
    with open(filename) as file_obj:
        for line in file_obj:
            if len(line) > 1:
                converted = list(line[:-1])
                # print(converted)
                list_of_list.append(converted)
    return list_of_list



