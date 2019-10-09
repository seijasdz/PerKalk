def divider(file):
    total = 0
    taas = 0
    tgas = 0
    tags = 0
    taa_matrix = []
    tga_matrix = []
    tag_matrix = []
    with open(file) as file:
        for line in file:
            if len(line) > 1:
                total += 1
                if line[8:11] == 'tag':
                    tags +=1
                    tag_matrix.append(list(line[:-1]))
                elif line[8:11] == 'tga':
                    tgas += 1
                    tga_matrix.append(list(line[:-1]))
                elif line[8:11] == 'taa':
                    taas += 1
                    taa_matrix.append(list(line[:-1]))

    taas_percent = taas / total
    tgas_percent = tgas / total
    tags_percent = tags / total
    print (taas_percent, tgas_percent, tags_percent)
    return taa_matrix, tga_matrix, tag_matrix

if __name__ == '__main__':
    divider('new_tts.exa')
