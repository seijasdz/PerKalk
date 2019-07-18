def get_codons(line):
    length = len(line)
    if 'atg' == line[0:3] and not length % 3:
        print('yepi')
        for idx, char in enumerate(line):
            if not idx % 3 and idx + 2 < length:
                pass
                print(char + line[idx + 1] + line[idx + 2])

    elif line[-3:] == 'tga':
        print('juejua')


def exon_trans_calc(filename):
    with open(filename) as file_handle:
        for line in file_handle:
            curated = line.replace(',', '').replace('p', '').replace('[', '').replace(']', '').replace('\n', '')
            print()
            get_codons(curated)

exon_trans_calc('exonesW.txt')
