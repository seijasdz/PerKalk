
def converter_to_two(sequence):
    new_list = []
    for index, element in enumerate(sequence):
        if index > 0:
            new_list.append(element + '|' + sequence[index - 1])
    return new_list

