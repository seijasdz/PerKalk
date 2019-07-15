import numpy


def create_matrix(filename):
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


def load_to_fit(filename):
    lines = []
    with open(filename) as file_handle:
        for line in file_handle:
            if line[0] != '>':
                lines.append(line[0:-1])

    return lines


lines = load_to_fit('TSS.f')
print(lines)