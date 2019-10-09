def seqs_from(filename):
    subseqs = []
    with open(filename) as file_handle:
        og_seq = ''
        for line in file_handle:
            tokenized = line.split()
            og_seq += ''.join(tokenized[0:-1])

        for string in (og_seq.split('N')):
            if string:
                subseqs.append(string.lower())
    return subseqs
