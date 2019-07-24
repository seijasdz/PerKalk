def to_string(file):
    with file.open() as file_handle:
        string = ''
        for line in file_handle:
            string += ''.join(line.split()[:-1])
        return string


def to_string2(file):
    with open(file) as file_handle:
        string = ''
        for line in file_handle:
            string += ''.join(line.split()[:-1])
        return string.lower()

