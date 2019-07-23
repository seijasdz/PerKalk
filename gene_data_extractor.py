from pathlib import Path
from xml.etree import ElementTree

genes = {}


def to_string(file):
    with file.open() as file_handle:
        string = ''
        for line in file_handle:
            string += ''.join(line.split()[:-1])
        return string


def slicer(string, divisions, complement, before=25, after=25):
    def reverse(s):
        return s[::-1]

    def first_replace(s):
        s = s.replace('A', 'W').replace('C', 'X').replace('G', 'Y').replace('T', 'Z')
        return s

    def second_replace(s):
        s = s.replace('W', 'T').replace('X', 'G').replace('Y', 'C').replace('Z', 'A')
        return s

    def correct(parts):
        corrected = []
        for part in parts:
            corrected.append(second_replace(first_replace(reverse(part))))
        return corrected

    data = {
        'complement': complement,
        'pure_cuts': [],
        'before_zone': [],
        'after_zone': [],
        'introns': []
    }

    for i, div in enumerate(divisions):
        data['pure_cuts'].append(string[div[0]:div[1]])

        before_cut = div[0] - before if div[0] - before > 0 else 0
        before_string = string[before_cut:div[0]] + 'P' + string[div[0]:div[1]]
        data['before_zone'].append(before_string)

        after_cut = div[1] + after if div[1] + after < len(string) else len(string) - 1
        after_string = string[div[0]:div[1]] + 'P' + string[div[1]:after_cut]
        data['after_zone'].append(after_string)

        if i < len(divisions) - 1:
            if div[1] < divisions[i + 1][0] and len(string[div[1]: divisions[i + 1][0]]) > 2:
                data['introns'].append(string[div[1]: divisions[i + 1][0]])

    if complement:
        data['pure_cuts'] = correct(data['pure_cuts'])
        c_after = correct(data['before_zone'])
        c_before = correct(data['after_zone'])
        data['after_zone'] = c_after
        data['before_zone'] = c_before

    return data


def transform(section):
    tokens = section.split('..')
    real_section = (int(tokens[0]) - 1, int(tokens[1]))
    return real_section


def get_divisions(sections):
    cut_pairs = []
    complement = False
    if 'A' in sections or 'Z' in sections:
        pass
    elif 'complement' in sections:
        print('complement')
        sections = sections.replace('join', '').replace('(', '').replace(')', '').replace('complement', '')
        pre_tuples = sections.split(',')
        complement = True
        for pre_tuple in pre_tuples:
            real_tuple = transform(pre_tuple)
            cut_pairs.append(real_tuple)
    elif 'join' not in sections:
        real_tuple = transform(sections)
        cut_pairs.append(real_tuple)
    else:
        sections = sections.replace('join', '').replace('(', '').replace(')', '')
        pre_tuples = sections.split(',')
        for pre_tuple in pre_tuples:
            real_tuple = transform(pre_tuple)
            cut_pairs.append(real_tuple)
    return complement, cut_pairs


def validate_cds(parts):
    if parts and parts[0][:3] != 'ATG':
        print('invalid start')
        return False
    if parts and parts[-1][-3:] != 'TGA' and parts[-1][-3:] != 'TAA' and parts[-1][-3:] != 'TAG':
        print('invalid end')
        return False
    return True


def divider(elements, sequence, tag='mRNA'):
    local = {}
    for el in elements:
        tokens = el.text.replace('\n', '').replace('\t', '').split('/')
        local[tokens[1]] = 'v'
        if tokens[1] not in genes:
            complement, divisions = get_divisions(tokens[0])
            data = slicer(sequence, divisions, complement)
            if tag == 'CDS':
                if validate_cds(data['pure_cuts']) and len(data['pure_cuts']):
                    print('len', len(data['pure_cuts']))
                    genes[tokens[1]] = data
            else:
                genes[tokens[1]] = data

    print(local)


folder_path = '/run/media/jose/BE96A68C96A6452D/Asi/Data/'
folder_path = '/run/media/zippyttech/BE96A68C96A6452D/Asi/Data/'

lookfor = 'CDS'
path = Path(folder_path)

subfolders = [x for x in path.iterdir() if x.is_dir()]

for folder in subfolders:
    xml = folder / 'Header.ebi.xml'
    sequence_body = folder / 'sequence_body.ebi'
    if xml.exists() and sequence_body.exists():
        print(folder)
        file_string = to_string(sequence_body)
        with xml.open() as file:
            patch = '<ROOT>' + file.read() + '</ROOT>'
        root = ElementTree.fromstring(patch)
        for ft in root.findall('FT'):
            elements = ft.findall(lookfor)
            divider(elements, file_string, lookfor)

pure_cuts = []
before_cuts = []
after_cuts = []
introns = []

pure_cuts_c = []
before_cuts_c = []
after_cuts_c = []


for key, gene in genes.items():
    if not gene['complement']:
        pure_cuts.append(' '.join(gene['pure_cuts']) + '\n')
        before_cuts.append(' '.join(gene['before_zone']) + '\n')
        after_cuts.append(' '.join(gene['after_zone']) + '\n')
        introns.append(' '.join(gene['introns']) + '\n')
    else:
        pure_cuts_c.append(' '.join(gene['pure_cuts']) + '\n')
        before_cuts_c.append(' '.join(gene['before_zone']) + '\n')
        after_cuts_c.append(' '.join(gene['after_zone']) + '\n')


pure_cuts = list(set(pure_cuts))
before_cuts = list(set(before_cuts))
after_cuts = list(set(after_cuts))

pure_cuts_c = list(set(pure_cuts_c))
before_cuts_c = list(set(before_cuts_c))
after_cuts_c = list(set(after_cuts_c))


with open('cuts.txt', 'w') as cuts, open('cutsa.txt', 'w') as acuts, open('cutsb.txt', 'w') as bcuts, open('cuts_intron.txt', 'w') as intron_file:
    for cut in pure_cuts:
        cuts.write(cut)

    for bcut in before_cuts:
        bcuts.write(bcut)

    for acut in after_cuts:
        acuts.write(acut)

    for gene_intron in introns:
        singles = gene_intron.split(' ')
        for c, i in enumerate(singles):
            if c > 2:
                if '\n' not in i:
                    intron_file.write(i + '\n')
                else:
                    intron_file.write(i)


with open('ccuts.txt', 'w') as cuts, open('ccutsa.txt', 'w') as acuts, open('ccutsb.txt', 'w') as bcuts:
    for cut in pure_cuts_c:
        cuts.write(cut)

    for bcut in before_cuts_c:
        bcuts.write(bcut)

    for acut in after_cuts_c:
        acuts.write(acut)
