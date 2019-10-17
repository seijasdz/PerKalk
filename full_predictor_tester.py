from pathlib import Path
from xml.etree import ElementTree
from gene_ebi_to_string import to_string


def reverse(s):
    return s[::-1]


def first_replace(s):
    s = s.replace('A', 'W').replace('C', 'X').replace('G', 'Y').replace('T', 'Z')
    return s


def second_replace(s):
    s = s.replace('W', 'T').replace('X', 'G').replace('Y', 'C').replace('Z', 'A')
    return s


def get_complete_cuts(divs, before, after, complement=False):
    new_divs = []
    for i, div in enumerate(divs):
        if i > 0:
            if not complement:
                intron = (divs[i - 1][1], div[0])
            else:
                intron = (div[1], divs[i - 1][0])
            new_divs.append((intron, 'n'))
        new_divs.append((div, 'an'))
    new_divs.insert(0, ((divs[0][0] - before, divs[0][0]), 'b'))
    new_divs.append(((divs[-1][1], divs[-1][1] + after), 'a'))
    return new_divs


def classify_characters(complete_cuts, string, complement=False):
    parallel = {
        'A': 'T',
        'C': 'G',
        'G': 'C',
        'T': 'A'
    }

    classified = []
    for cut in complete_cuts:
        start = cut[0][0]
        end = cut[0][1]
        if not complement:
            for base in string[start:end]:
                classified.append((base, cut[1]))
        else:
            for base in string[start:end][::-1]:
                classified.append((parallel[base], cut[1]))

    if complement:
        pass
        # print(complete_cuts)
    return classified


def slicer(string, divisions, complement, before=25, after=25):
    if not complement:
        if divisions:
            classified_chars = classify_characters(get_complete_cuts(divisions, before, after), string)
            start = divisions[0][0] - before
            end = divisions[-1][1] + after
            section_with_gene = string[start:end]

            pure_cuts = [string[div[0]: div[1]] for div in divisions]
    else:
        if divisions:
            if len(divisions) > 1:
                start = divisions[-1][1] - after
                end = divisions[0][0] + before
            else:
                start = divisions[0][0] - after
                end = divisions[-1][1] + before

            section_with_gene = second_replace(first_replace(reverse(string[start:end])))
            classified_chars = classify_characters(get_complete_cuts(divisions, after, before, complement), string, True)
            print(section_with_gene)
            print(divisions[:10])
            print(classified_chars[:28])

def transform(section):
    tokens = section.split('..')
    real_section = (int(tokens[0]) - 1, int(tokens[1]))
    return real_section


def get_divisions(cut_annotations):
    cut_pairs = []
    complement = False
    if any(x in cut_annotations for x in ['A', 'Z', 'U']):
        return complement, cut_pairs
    elif 'complement' in cut_annotations:
        complement = True
    cut_annotations = cut_annotations.replace('join', '').replace('(', '').replace(')', '').replace('complement', '')
    pre_tuples = cut_annotations.split(',')
    for pre_tuple in pre_tuples:
        real_tuple = transform(pre_tuple)
        cut_pairs.append(real_tuple)
    return complement, cut_pairs


def cut(element, file_string):
    tokens = element.text.replace('\n', '').replace('\t', '').split('/')
    cut_annotations = tokens[0]
    gene_id = tokens[1]
    complement, cut_pairs = get_divisions(cut_annotations)
    slicer(file_string, cut_pairs, complement)


def get_lines(annotations, file_string):
    return [cut(annotation, file_string) for annotation in annotations]


def test(tag, folder):
    path = Path(folder)
    sub_folders = [x for x in path.iterdir() if x.is_dir()]
    datas = []
    for folder in sub_folders:
        xml = folder / 'Header.ebi.xml'
        sequence_body = folder / 'sequence_body.ebi'
        if xml.exists() and sequence_body.exists():
            file_string = to_string(sequence_body)
            with xml.open() as file:
                patch = '<ROOT>' + file.read() + '</ROOT>'
            root = ElementTree.fromstring(patch)
            for ft in root.findall('FT'):
                annotations = ft.findall(tag)
                datas.append(
                    {
                        'annotations': annotations,
                        'file_string': file_string
                    }
                )
    lines = [get_lines(x['annotations'], x['file_string']) for x in datas]

if __name__ == '__main__':
    route = '/run/media/jose/BE96A68C96A6452D/Asi/Data/'
    # route = '/run/media/zippyttech/BE96A68C96A6452D/Asi/Data/'
    test('CDS', route)
