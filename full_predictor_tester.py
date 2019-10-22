from pathlib import Path
from xml.etree import ElementTree
from gene_ebi_to_string import to_string
from converter_to import converter_to
import time
import numpy
from gene_predictor import predict_all

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
        new_divs.append((div, 'a'))
    if not complement:
        new_divs.insert(0, ((divs[0][0] - before, divs[0][0]), 'b'))
        new_divs.append(((divs[-1][1], divs[-1][1] + after), 'f'))
    else:
        new_divs.insert(0, ((divs[0][1], divs[0][1] + before), 'b'))
        new_divs.append(((divs[-1][0] - after, divs[-1][0]), 'f'))
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
    return classified


def slicer(string, divisions, complement, before=25, after=25):
    section_with_gene = None
    classified_chars = None
    if not complement:
        if divisions:
            start = divisions[0][0] - before
            end = divisions[-1][1] + after
            new_divs = get_complete_cuts(divisions, before, after)
            section_with_gene = string[start:end]
            classified_chars = classify_characters(new_divs, string)
            # pure_cuts = [string[div[0]: div[1]] for div in divisions]

    else:
        if divisions:
            start = divisions[-1][0] - after
            end = divisions[0][1] + before
            new_divs = get_complete_cuts(divisions, after, before, complement)
            section_with_gene = second_replace(first_replace(reverse(string[start:end])))
            classified_chars = classify_characters(new_divs, string, True)

    return {
        'section_with_gene': section_with_gene,
        'classified_chars': classified_chars
    }


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
    sliced = slicer(file_string, cut_pairs, complement)
    if sliced['classified_chars']:
        sss = [''.join(x) for x in sliced['classified_chars']]
        return ' '.join(sss) + '-' + str(complement)

def get_lines(annotations, file_string):
    return [cut(annotation, file_string) for annotation in annotations]


def get_annotated_lines(tag, folder):
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

    return lines


def get_testable_string(annotated_string):
    l = [x[0] for x in annotated_string.lower().split()]
    a = [x[1] for c, x in enumerate(annotated_string.split()) if c > 1]
    two = converter_to(l, 2)
    seq = numpy.array(two, numpy.unicode_)
    return seq, a,''.join(l)



if __name__ == '__main__':
    route = '/run/media/jose/BE96A68C96A6452D/Asi/Data/'
    # route = '/run/media/zippyttech/BE96A68C96A6452D/Asi/Data/'

    anno = "Cb Ab Cb Cb Tb Cb Tb Tb Cb Tb Gb Cb Cb Ab Cb Ab Ab Ab Cb Gb Tb Cb Ab Gb Cb Aa Ta Ga"
    seq, a, string = get_testable_string(anno)
    print(seq)
    print(a)
    print(string)

    predict_all(seq, string)

    #lines = get_annotated_lines('CDS', route)

    #with open('out_CDS_annotated.txt', 'w') as o_file:
    #    for l in lines:
    #        for line in l:
    #            if line:
    #                o_file.write(line + '\n')
