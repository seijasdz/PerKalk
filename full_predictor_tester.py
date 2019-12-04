from pathlib import Path
from xml.etree import ElementTree
from gene_ebi_to_string import to_string
from converter_to import converter_to
import time
import numpy
from gene_predictor import predict_all
from math import sqrt

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


def cut(element, file_string, part):
    tokens = element.text.replace('\n', '').replace('\t', '').split('/')
    cut_annotations = tokens[0]
    gene_id = tokens[1]
    complement, cut_pairs = get_divisions(cut_annotations)
    sliced = slicer(file_string, cut_pairs, complement)
    if part == 'complement' and complement:
        if sliced['classified_chars']:
            sss = [''.join(x) for x in sliced['classified_chars']]
            return ' '.join(sss) + '-' + gene_id
    #else:
    #    if sliced['classified_chars']:
    #        sss = [''.join(x) for x in sliced['classified_chars']]
    #        return ' '.join(sss) + '-' + gene_id

def get_lines(annotations, file_string, part):
    return [cut(annotation, file_string, part) for annotation in annotations]


def get_annotated_lines(tag, folder, part):
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
    lines = [get_lines(x['annotations'], x['file_string'], part) for x in datas]

    return lines


def get_testable_string(annotated_string):
    l = [x[0] for x in annotated_string.lower().split()]
    a = [x[1] for c, x in enumerate(annotated_string.split()) if c > 1]
    two = converter_to(l, 2)
    seq = numpy.array(two, numpy.unicode_)
    return seq, a,''.join(l)


def calculate_accuracy(ann, path):
    if len(ann) != len(path) - 1:
        print("incompatibles lengths!!")
        return

    an = ['start zone8', 'start zone9', 'start zone10', 'start zone11',
          'start zone12', 'start zone13', 'start zone14', 'start zone15',
          'start zone16', 'coding state 0', 'coding state 1', 'coding state 2',
          'donor10', 'donor11', 'donor12', 'donor13',
          'acceptor116', 'acceptor117', 'acceptor118', 'acceptor119', 'acceptor120',
          'donor00', 'donor01', 'donor02', 'acceptor016', 'acceptor017', 'acceptor018',
          'donor20', 'donor21', 'donor22', 'donor23' 'donor24',
          'acceptor216', 'acceptor217', 'acceptor218', 'acceptor219',
          'stop zone taa8', 'stop zone taa7', 'stop zone taa6', 'stop zone taa5', 'stop zone taa4',
          'stop zone taa3', 'stop zone taa2', 'stop zone taa1', 'stop zone taa0',
          'stop zone tag8', 'stop zone tag7', 'stop zone tag6', 'stop zone tag5', 'stop zone tag4',
          'stop zone tag3', 'stop zone tag2', 'stop zone tag1', 'stop zone tag0',
          'stop zone tga8', 'stop zone tga7', 'stop zone tga6', 'stop zone tga5', 'stop zone tga4',
          'stop zone tga3', 'stop zone tga2', 'stop zone tga1', 'stop zone tga0']

    def zone_counter():
        return {
            'tp': 0,
            'tn': 0,
            'fp': 0,
            'fn': 0
        }

    useful_path = [c for i, c in enumerate(path) if i]
    true_positives = 0
    fake_positives = 0
    true_negatives = 0
    fake_negatives = 0
    zone_counters = {
        'start': zone_counter(),
        'coding': zone_counter(),
        'donor': zone_counter(),
        'acceptor': zone_counter(),
        'stop': zone_counter(),
    }

    """-start zone

start zone12
start zone13
start zone14

-donor

donor14
donor15
donor03
donor04
donor25
donor26

-accep

acceptor121
acceptor122
acceptor019
acceptor020
acceptor220
acceptor221

stop zone"""

    for i, c in enumerate(ann):
        if c == 'b' or c == 'n' or c == 'f':
            if useful_path[i] not in an:
                true_negatives += 1
            else:
                fake_positives += 1
        else:
            if useful_path[i] in an:
                true_positives += 1
            else:
                fake_negatives += 1

    try:
        sensitivity = true_positives / (true_positives + fake_negatives)
        specificity = true_positives / (true_positives + fake_positives)
        correlation_coefficient = ((true_positives * true_negatives) - (fake_negatives * fake_positives)) / \
            sqrt((true_positives + fake_negatives) * (true_negatives + fake_positives)
                 * (true_positives + fake_positives) * (true_negatives + fake_negatives))
        return {
            'sensitivity': sensitivity,
            'specificity': specificity,
            'cc': correlation_coefficient,
            'true_positives': true_positives,
            'true_negatives': true_negatives,
            'fake_positives': fake_positives,
            'fake_negatives': fake_negatives
        }
    except ZeroDivisionError:
        print('--invalid-- zero division cc')
if __name__ == '__main__':
    # route = '/run/media/jose/BE96A68C96A6452D/Asi/Data/'
    route = '/run/media/zippyttech/BE96A68C96A6452D/Asi/Data/'

    #anno = "Cb Tb Gb Cb Tb Cb Cb Cb Ab Cb Cb Ab Gb Cb Tb Tb Gb Gb Ab Gb Ab Cb Ab Ab Cb Aa Ta Ga Ta Ga Ga Ta Ta Ca Ta Ta Ga Aa Ca Aa Aa Ca Ta Ca Ta Ga Ca Ta Ca Ca Ta Ta Ta Ga Ga Ga Gn Tn An An Gn Tn Tn Gn Gn An Cn Tn Cn An Gn An Gn Gn Gn Gn An Cn An Gn Tn Tn An Gn An An Gn Gn Gn Tn An Cn An Gn Gn Cn Tn Gn Tn Gn Gn Cn Tn Gn Tn Tn Gn Tn Gn An Gn Tn Cn An An Gn An Gn Tn Tn Tn Tn Gn Tn Cn Tn Tn Cn Cn Tn Gn Tn Gn Gn Tn An An Cn Tn Cn Tn Gn Gn Gn Tn An Gn An An Cn Tn Cn An Tn Gn An Gn Tn An Tn Gn An An Gn Cn An An Cn Tn Tn Gn Tn An Tn Cn Tn Gn Tn Gn Cn Tn Tn Cn Cn An Tn Gn Gn Tn Tn Tn An Tn Tn An Gn An Gn Cn Tn Tn An Tn Tn Tn Tn An Tn Gn An An An An Gn Gn An Tn Gn Gn Gn An An Gn Gn Gn Cn An An Cn Cn Cn Tn Gn An Gn Gn Tn An Gn Cn An Tn Tn An An Gn Cn Cn Tn Gn Gn An Cn Gn Cn An Cn Cn Gn Cn An Gn Tn Gn An An Gn Tn Tn Tn Cn Cn Tn Tn Gn An Tn An An Cn Cn An Cn Cn Tn Gn Tn An Gn Cn Tn Tn Gn Tn Tn Cn An Gn Tn Tn Cn Tn Gn Tn Tn An Gn Tn An Cn Tn Gn Gn An Tn Tn Tn Tn Gn An Gn An An An Gn An Gn An An An Tn An Gn An An An Cn Tn Cn An An Gn An Gn An Tn Cn Tn Gn An Gn Tn Tn Gn An Tn Cn Cn Cn Tn Cn An Gn An Gn Tn Cn Tn An Cn An Tn Tn An An Tn Tn Cn Tn Gn Tn Cn Tn Cn Cn Cn Cn An An Tn Tn Cn Tn Cn Tn Cn Tn Tn Cn Cn Tn Cn An Tn Tn An Tn Tn Tn Tn Cn Cn Tn Tn Gn Gn An Cn Cn An An Cn Tn Gn An Tn An Tn Cn Tn Tn Tn An Tn Tn Cn Tn Cn Tn Gn An Tn Cn Tn Cn Tn Tn Gn Cn An Gn Tn Tn Cn Cn An Gn Tn Tn Gn An Tn Gn Gn Gn Cn An An Gn Tn Gn Gn Gn Tn Gn An Gn Tn Gn An Tn Cn Tn Cn Tn An An Cn Tn Cn An Gn Cn Tn Tn Cn Tn Cn Cn Tn Tn Cn Tn An Tn Gn Cn Cn An Cn Tn Tn Tn Cn Cn Tn An Cn Tn Tn Cn Cn An An An Gn Gn An Tn Gn Gn Gn Tn Cn Cn Tn An Tn Tn An An Cn Cn Tn Gn Cn An Gn An An Gn An Gn Cn An Tn An Tn An Gn Gn Gn An An An Gn Cn An Gn An Gn An An An Gn An An Gn An An An Gn An Tn Tn Tn An Tn An An An Tn Tn An Tn Gn Tn An An An An Tn Cn Tn Cn Cn Cn An Tn An Tn Tn Cn An Gn An Gn Cn An Tn Gn An Tn Tn Tn Tn Cn An Tn An An An An An Cn Tn An Tn Gn An Tn Tn Tn Gn Tn An An Gn Cn Tn Tn Tn Cn Tn Tn Tn Cn Tn An An An Gn Tn Tn Tn Tn Cn Tn Gn An An Tn Gn An Gn Cn Tn Cn Tn Tn Tn Cn Cn Tn Cn Tn Gn An Gn Cn Gn An An An Tn Gn Cn Cn Tn Gn An An Tn Gn Tn An Tn An Tn An Tn An An Cn An Gn Cn An Gn Gn An An An Gn Gn Gn Cn Cn An Gn An An Tn Tn Gn Cn An Cn An An An Cn An Cn Tn An An An Gn Tn An Gn Tn Cn Tn Gn Cn Tn Gn An Cn An Cn Tn Gn An An Gn Tn An Gn Tn Cn Tn Gn Cn Tn Gn Gn Gn Tn Tn Tn Cn Cn Cn An An Cn An Tn Cn Cn Cn Cn Tn An An An Cn Cn An Gn Tn Cn Cn Cn Cn Tn Tn Tn An Cn An An Gn An An An Gn Cn An Gn An Gn Tn Gn Gn Gn Tn Tn Tn Tn Tn Tn Cn Tn Cn Cn Tn Cn An Gn An An An Gn An Gn Cn Cn An Gn Gn Cn An An An Tn An Cn Tn Gn An An Gn Gn Cn Tn Tn Cn An An An Gn An Gn An Gn Gn An Gn Tn Tn Tn Gn An An Cn Cn Tn Cn Cn Gn Tn Tn Tn Cn Cn An An Tn Gn Tn An Cn An An Gn Cn Tn Gn Cn Tn Gn An An Gn Cn Cn Cn Tn Gn Cn An An Tn Gn Tn An Cn Tn Cn Tn Cn An Gn An Gn An Cn An Gn Tn An Tn An Gn Cn Tn Cn An Gn An Gn Gn An Tn An An Gn Gn An An An An Cn An Gn Cn Tn Cn Cn Gn Cn Tn Tn Gn Gn Cn Tn Tn Gn Cn Cn Cn An An Gn Tn Cn Tn An An An An Cn Cn Cn Tn Gn Cn Tn Tn Cn An Tn Tn Gn Cn Tn Cn Cn Cn An An Cn An Gn An An An An An Cn Cn Tn An Cn An Tn Cn An Gn Tn An An Tn Tn Tn An An Tn Cn Cn Tn Cn An Gn Gn An Gn Cn Cn Cn Cn An Gn An Tn Tn Cn Cn Cn An An Cn An Cn Cn Cn An Gn An Tn Gn An Gn Gn An Tn Tn Tn An An Tn An An Cn An An Tn Gn An An Cn An Gn An Gn Gn Tn Cn Tn Cn Tn Gn Cn An Gn Cn An Gn Gn An Gn Gn An Cn An Cn Tn Cn Tn Tn An An An Gn Gn Cn Tn Gn Gn Cn Tn Tn Tn Tn Gn Tn An Gn Gn Tn Tn Tn Cn Cn An Cn Tn Tn Gn Cn Cn Tn An Tn Tn Tn Tn Tn Cn An An An Tn An An An Tn Cn Tn An An An Tn An Tn Gn Tn An Tn Tn Tn Tn Tn Cn Cn Cn Tn Tn Gn Tn Gn Cn Tn An Tn Cn Cn An Cn An An An An Gn Cn Cn Tn Cn An Cn Cn An Gn Tn Tn Gn Cn Tn Cn An Gn An Gn Cn Cn Cn Tn Cn Cn Cn Cn Tn Tn Tn An Cn Tn Cn Cn An Cn Tn Tn Gn An Tn Cn Tn Tn Tn Tn Tn Cn Tn Tn Tn Tn Tn Cn Tn Tn Cn Tn Tn Gn Tn An An Tn Cn Tn Cn Cn An An Gn Tn An Gn An Cn An Cn Cn An Cn An An An Gn Gn Cn An Gn Tn Gn An Tn Cn An Cn Tn Tn Tn Gn Cn An Gn Cn Cn Tn Cn Cn An Tn Gn Gn Gn Tn Cn An Gn Cn Gn Tn Gn Tn Tn Cn Cn An An Gn An Gn Gn An An An Cn Cn Gn Tn An An Cn Cn Tn Tn Gn Cn An Cn Tn Gn Tn Gn An Gn Gn Tn Gn Cn Tn Cn Cn An Tn Cn Tn Gn Cn Cn Tn Gn Gn Gn An Gn Cn An Gn Cn Tn Cn Tn An Cn An Cn An Gn Tn Gn Gn Tn Tn Tn Cn Tn Cn An An Tn Gn Gn Cn An Cn An Gn Cn Cn An Cn Tn Cn An Gn An Cn Cn Tn Cn Gn An Cn Cn Cn Cn Cn An Gn Cn Tn An Cn An Gn An An Tn Cn An Cn Cn Tn Cn Tn Gn Cn Cn An Gn Tn Gn Tn Cn An An Tn Gn An Cn An Gn Tn Gn Gn Tn Gn An An Tn An Cn An Gn Gn Tn Gn Cn Cn An Gn An Gn An Gn Gn Tn Cn Tn Cn Tn Cn An Gn Gn Gn Cn Gn An An Gn Tn Gn An Cn Cn Cn Cn An Tn An Cn An Gn Cn Tn Gn Gn An An An Tn Cn Cn An Cn An Gn An Gn Gn Tn An An Tn Tn An Tn Gn An Cn Tn Tn Gn Gn An Cn Cn An Gn Gn An Gn Gn Gn Cn Cn Gn Gn An An An Cn Cn An Cn An An Gn Gn Tn Cn Tn Tn Cn Cn Tn Cn Tn Gn Cn Gn Tn Tn Cn Tn Cn Tn Cn Cn Tn Tn Tn Cn Tn Gn Gn An Tn Gn Cn Cn An Gn An Tn Gn Tn Gn Tn Gn Tn Gn Tn Gn Tn Cn Gn An Tn Tn Tn An Tn Cn Tn Gn Gn Gn Tn Gn Cn An Gn Gn Tn Cn Tn Cn An Gn Cn An Cn Tn An Tn Gn Tn An Cn Cn Tn An Cn Cn An Gn Gn Tn Gn An An An Cn Tn An An An Tn Cn An Gn Tn An Tn An Cn Tn Cn An An An An Cn An Tn Cn An Cn An Gn Cn An Gn Gn Gn Tn Cn Cn Tn Tn An An An Cn Tn Gn Tn Gn Tn An Cn Cn An An Cn An Gn Tn An Cn Cn Tn Cn An Gn Cn Tn An Gn Gn Cn Cn An Gn Cn Tn Cn Tn An An Gn Tn An Gn Gn An Tn Gn Tn Cn Tn An Gn An Gn Tn An An Gn An Cn Tn Tn An Gn Gn Cn Cn Cn Tn An An Gn An Tn Tn Gn Tn An Tn Cn An Tn Tn An Tn Tn Cn Cn An An An Tn Gn Tn An Tn Gn Cn Cn Tn Gn Tn An Cn An Tn An Gn Tn An An An Tn An An Tn An An Tn An An Tn Tn Cn An An An Cn An An An Tn Gn An Gn Cn An An Tn An An Tn Gn An An Cn An An Tn Gn An An An Tn An Tn Cn Tn Tn Cn Cn Tn Cn Cn Cn An Cn Cn Cn An An An An Cn Cn Tn An Gn An Tn Cn Tn Gn Cn An Gn Tn Cn An Cn Cn Tn Tn Cn Cn An Cn An Gn An An An Cn An An An An Cn Tn Gn Cn Tn An An Tn An Tn An Tn An Tn An Tn An Cn An Tn An Tn Tn An Tn Gn Cn An Tn An Tn An Tn An Tn Gn Tn An Tn An Tn An Tn Tn An Gn Cn An Tn Tn An Tn Tn An Tn An Tn An Tn An An An Cn Tn An Tn An Tn An Tn An Tn An Tn An Tn An Tn An Tn An Tn An Cn An Cn An Cn An Tn An Tn An Tn An Tn An Gn Tn Tn Tn Tn An An An Gn Cn Tn Gn Tn Gn Gn Gn Tn Gn Gn Tn An Gn An An Cn An Cn Tn Tn Tn Cn Tn An An Gn Tn An Tn Tn Gn Tn Cn Cn Tn An Tn An Tn Cn Tn Tn An Cn Tn Tn Tn Tn An Tn Cn An Cn Tn Tn An An An An Tn An Cn Cn Tn Tn An An An Gn Tn Tn Cn Tn Gn Tn Tn Tn An Tn An Tn Tn An An Tn An Gn An An An Tn Cn Cn An Cn Tn Tn Tn An An An An Tn Gn Tn Cn Tn Gn Cn Tn Gn Cn An Cn An Gn Cn Tn Gn Tn An Cn Cn An Gn An An Tn Tn An An Tn Tn Tn An An Cn Cn An Tn Tn Cn Tn Gn Cn Cn Tn An Tn Tn An Gn Cn An An An Tn Tn Tn An Tn An An Gn An Tn Gn Tn Tn Tn Cn Cn An Cn Tn Cn Tn Tn Tn Tn Gn An Tn Tn Tn Tn An Cn An An Gn Cn An Tn Tn Gn Cn Tn Gn Cn An An Tn An An Tn Tn An Tn Cn Tn Tn Tn An Tn An Cn An Tn An Cn An Tn Tn Cn Cn Gn Gn Gn Gn Gn Cn An Cn An Tn Tn Tn Cn An Tn An Gn Gn Gn Tn An Tn Cn Tn Gn Tn An An Gn An Tn An An Tn An Tn Tn Cn Tn Gn Gn An Cn An Tn Gn Gn An An Tn Tn An Tn Tn Gn Tn Cn An Tn An Tn Tn Tn Tn An An Tn An An An Tn Gn Gn An An Gn Cn An Gn An An An An Gn Cn Cn Cn An Cn Cn An Cn Cn Cn An Cn Tn Tn An An An Cn Cn Tn Gn Tn An An Gn Gn An Gn Tn Tn Tn Cn Tn An Tn Tn Cn An Tn Cn An An Gn An Gn An Tn Gn Gn An Gn An Cn Tn Gn Cn An An Tn An Tn Tn Cn Tn Tn Cn An Gn Gn Gn An An An An An Tn Gn Cn Tn Tn Cn An Tn Tn Cn An An Tn An Tn Gn Cn An An Gn An Cn An An An An Tn An An Cn Cn An An An Tn Gn Tn An Tn Tn An Tn Tn Tn An Gn Tn Tn An Tn Tn Cn An Gn An An Gn Tn Gn An Tn Gn An Gn Tn Tn Tn An An An Cn Tn Tn Cn Gn Tn Tn An An An An An Tn Gn An Tn An An Tn Gn Cn An An An An Gn Tn Cn An An Tn Cn Cn An Cn An Cn Tn Gn Gn An Gn Gn Cn An An An An Tn Tn Gn Cn Cn Tn Tn Tn An Tn Gn Cn An Gn An Cn An Tn Gn Gn Tn Tn Gn Cn Tn Gn Tn Tn An Gn An Gn Gn An An An Tn Gn Tn Gn Tn Tn Tn Tn An Gn Cn Cn An An Gn Gn An An An Tn Gn Tn Tn Gn An Gn An An Tn Gn Tn Cn Tn Gn An Gn Gn Gn An An Gn Gn Cn Tn Tn Gn An Tn Gn Tn Tn Cn Cn Tn Cn An Tn Gn Tn Gn Tn Cn Cn Cn Tn Tn An Gn Gn Tn Cn Tn Cn Tn Gn Tn An An An Tn An An Gn Tn Tn Cn Tn Tn An Gn Gn An Tn Gn Gn An Gn An Gn Gn Tn Tn An Gn Gn Gn An An Gn Cn Tn Gn An Cn An Gn An Gn Cn Tn Gn Tn Tn Tn Cn Gn Tn Tn Tn Tn Tn Tn Tn Tn Tn Tn Tn Tn Tn Tn Tn Tn Tn Tn Tn Tn Tn Tn Tn Tn Tn Tn Gn An Gn An Cn An Gn An Gn Tn Cn Tn Cn An Cn Tn Cn Cn Gn Tn Tn Gn Cn Cn Cn An Gn Gn Cn Tn Gn Gn An An Tn Gn Cn An Gn Tn Gn Gn Tn Gn Tn Gn An Tn Cn Tn Tn Gn Gn Cn Tn Cn An Cn Tn Gn Cn An An Cn Cn Tn Cn Tn Gn An Cn An Gn An Gn Cn Tn Gn Tn Tn Tn Cn An An Cn An An Gn An Tn An Cn Tn Tn Tn Gn Gn Gn Tn An An Tn Tn Cn An An Gn Cn An Gn An An Gn Tn Gn Tn Tn Tn Cn An Tn Gn Gn Tn Tn Cn An An Gn Tn Gn An Cn An An Cn An An An An Gn Gn An Cn Tn An Gn Gn Gn Cn Tn Gn Gn Gn An An Gn An Gn Tn Cn Cn An Tn Tn Tn Tn Cn Tn Tn Tn An Cn Cn Tn Cn Tn An Cn An Cn An Gn An Cn Tn Cn Tn An Gn Gn An Cn An Cn Tn Cn An Gn Cn An Tn Cn An Gn Tn Gn An Cn Tn Tn An Cn Tn Gn An Gn An An Tn Cn An Tn An Tn Cn Tn Tn Tn Tn Tn Tn An Tn An Tn Tn An Cn An Tn Tn Tn Tn An An Tn Gn Tn Gn An An Tn Tn Tn Tn Cn An Tn Tn Gn An An An Tn An Tn An An Cn Tn Tn An Tn Gn Tn An Cn An An Cn An An Gn Cn Cn An Tn An An An Tn Gn Tn An Tn Gn Cn Cn An Cn An Gn Tn An An An Tn Tn Tn Tn Tn Tn Cn An An An Gn Tn An An An Cn An Gn Cn An Tn Cn Cn An Gn An Tn Cn Tn An Tn Tn Tn Tn An Gn An Cn An Tn Tn Tn Cn Cn An Gn Cn An Cn Cn Tn Cn An Gn An An An Cn Tn Cn Cn Tn Tn Cn Tn Tn Gn Cn Cn Cn An Tn Tn Tn Cn Tn An Gn Tn Cn Cn An Cn Cn Cn Cn Tn Gn An Cn An Cn Tn Gn An Cn Tn Tn Cn Cn An Gn Cn An Cn Cn An Tn An Gn Gn Tn Tn Tn Tn Tn Tn Gn Cn Cn Tn An Tn Tn Tn Tn Tn Gn An An Cn Tn Tn Cn An Tn An Tn An An An Gn Gn Gn Cn Tn An Cn Gn Gn Tn An Tn Gn Gn Gn An Cn An Gn Gn An Tn Gn Tn An Cn Tn Cn Tn Tn Tn Tn Gn Tn Gn Tn Cn Tn Gn Gn Tn Tn Tn An Tn Tn Tn Tn Tn Cn Tn Cn An An Cn An Tn Tn An Tn Gn An Tn Tn Gn Tn Gn An An Cn Tn Gn Tn An Tn An Cn An Tn Gn Tn Tn Gn Tn An An An Tn Tn Gn An An Gn Tn Tn Cn An Tn Tn Tn Gn Tn Tn Cn Tn Tn An Tn Tn An Cn Tn Gn Tn An Tn An Tn Tn Tn Cn An Tn An Tn Tn Tn Tn Cn Cn Cn An Tn Tn Cn Tn Gn Cn Tn Gn An Gn Gn Gn Gn Cn Tn Cn Gn Gn Gn Tn Tn Gn Cn Tn Tn Cn Tn Gn Gn Tn Tn Tn Tn Gn An Cn Cn An Cn Tn Tn An Cn An An An Tn An Gn Tn An An Tn Gn Cn Tn Gn Tn Gn Gn An Tn An Tn Tn Tn Tn Tn Gn Tn Tn Tn An Cn An Tn Gn Tn Tn Tn Cn An An An Tn An Tn Tn Tn Tn Cn Tn An Cn Cn An Tn Tn An Cn Tn An An Tn Gn Tn Tn Tn Gn An Cn An Gn Tn An Tn Cn Tn Gn An Tn Gn Gn An Gn Gn Tn An Tn Gn Tn Tn An Cn An An Tn Tn Tn Cn Cn Cn An Cn Tn Cn Tn Tn Tn Tn Gn Gn An Tn An Tn An Cn An An An Tn Gn Cn An An An Tn Tn An An Gn Cn An Tn Tn Tn Cn Cn Cn Tn An An Tn Gn An Cn Tn An Tn Gn An Tn An Tn Tn Gn An Cn Cn An Tn Tn Tn Tn Tn Cn An Cn Gn Tn Tn Tn An Tn Tn An Cn Cn Tn An Tn Tn Tn Gn Gn Gn Tn An Tn Cn Cn Tn Tn Tn Tn Tn Tn Gn Tn An An Gn Gn Tn Cn Tn Tn Tn Tn Tn Gn An Tn An An Tn Tn Tn Tn Tn Cn Tn An Tn Tn Gn Tn Cn Tn Gn Tn Cn Tn Tn Tn Tn Tn Cn Tn Tn An Cn Tn Gn An Cn Tn Tn An Cn An Gn Gn An Gn Tn Tn An Tn Tn Tn An Tn An Tn An Tn Tn Cn Tn Gn Gn An Tn An Tn An An Gn Tn Cn Cn Tn Tn Tn Tn Tn Cn An An Gn Tn Gn Tn Cn Tn Tn An An An Tn An Tn Cn Tn Tn Cn Tn Cn Cn Cn An Cn Tn Cn Tn Gn Tn Gn Gn Cn Tn Tn Tn Tn Cn Tn Tn Tn Tn Cn An Cn Tn An Tn An An Tn An Tn Tn An An An Gn Tn Tn Tn Tn Gn Tn Tn Tn Gn Tn Tn Tn Tn Cn An An Tn Gn An An Cn An An Gn An Cn Tn Tn Gn Tn Tn An An Tn Tn Tn Tn An Cn Tn An Cn An Gn Tn An Cn An An Tn An Cn An Cn Cn An An Cn An Tn Tn Tn Tn Tn Gn Gn Tn Tn Gn Gn Tn Gn Cn Cn Gn Tn Cn Tn An Tn Gn Cn Cn Cn Tn An Tn Tn Tn Tn An An An An Tn Gn Tn Tn Tn Tn Gn Cn Cn An Cn An Gn Tn Cn An Gn Cn Cn An Tn Gn Gn Gn Gn An An Tn Tn Cn An An Tn Tn Gn Cn An Gn Cn Tn Cn An Cn Tn Gn An Cn Tn Cn Tn Tn An Gn Tn Cn An Tn Tn Tn Cn An An Tn Tn Tn An An Tn Cn An Gn Tn An Tn An Tn Gn Cn Tn An An Gn An Cn Cn An An Gn An An An An Tn Gn Cn An Gn Cn Tn Gn Gn Tn Tn Gn Gn Tn Gn Gn An Cn Cn Cn Gn An Gn Cn Cn Cn Tn An An Tn Cn An Tn An Tn Gn Gn Cn Tn Gn An Cn An Cn Cn An An An Gn Tn Cn An Tn Gn An An Gn Gn Tn An Gn Gn Gn Gn An An Tn Gn An Tn Tn Gn Tn Tn Gn An Gn Tn An Gn Gn An An An Tn Tn An Gn Tn An An Tn Gn Tn Cn Tn Gn Cn Tn An An An An Gn Cn Tn Gn Cn An Tn Cn Tn Cn Gn Tn Gn An Gn Tn Tn Tn Tn Gn An Tn An Tn An Tn An An Tn An Cn Tn An Cn Cn An Tn Tn An Cn Cn An Tn Tn Tn An An Gn Tn Tn Tn An An An An An Cn Tn Cn An Tn An An An Tn Tn Tn Cn Cn An Tn Tn An Tn An Cn An Tn Tn An Tn Tn Tn Tn An Tn Tn Cn An Tn Gn An An Tn Tn An Cn Tn Tn An Gn An Tn Gn Tn Tn Tn Tn Tn Tn Gn Gn An Gn Gn Gn Tn An Gn Tn Tn Tn Tn Tn Tn Gn Cn Tn An Tn Tn Gn An Tn Tn Tn Cn Tn An An Tn An Tn Tn An Tn An Gn Cn An Cn Tn Tn Tn Gn Gn Tn Tn An Gn An An Tn Tn Tn Gn Tn Gn Tn Tn An An An Cn An Tn Gn Gn Tn An Tn Cn An Gn Tn Tn Cn An Tn Tn An Gn Tn An Tn Tn Tn Gn Tn Tn Gn An An An An Tn Tn Gn Cn Tn Tn Tn Gn Tn Gn Gn Cn Cn Tn An Gn Tn Tn Cn An Cn An Cn Tn Tn An An Tn Tn Tn Tn Tn Cn Tn An An Tn An Tn Tn Gn Cn An Tn Gn Tn Gn Tn Gn Tn Tn Tn Gn An An An An Gn Cn An Tn An Tn Gn An An Tn Tn Tn Tn Cn Tn Tn Tn Cn Tn An An An Gn Tn An Gn An Gn Cn An Tn Cn Cn An Cn An An An An Tn Cn An An Gn Cn Cn Tn An Tn Tn An An Tn Tn Gn Tn Gn Tn Cn An Tn Tn Cn An An An Tn An Tn Tn Tn Tn Cn Tn An Cn Cn An Tn Tn An Cn Tn An Gn Tn Gn Tn Tn Tn Gn An Cn An Gn Tn Tn Tn Cn Tn Gn An Tn Gn Gn An Gn An Tn Gn Tn Gn An Cn An Gn Cn An Cn Tn Tn Cn Cn Cn An Cn Tn Cn Tn Tn Tn Tn Gn Gn An Tn An Tn An Tn An An An Gn Cn An An An Gn An Tn Gn Gn Cn An Gn Tn Gn Tn Tn Tn Gn Gn Gn Cn An Gn Tn Gn Cn An Gn Tn Cn Tn Tn Tn Cn An Gn Cn An Tn Cn Tn Gn Cn Tn Gn Gn Cn Tn Gn Tn Gn Gn Gn Tn Gn Gn Gn An An Tn An An Gn Gn Gn An Gn Gn Gn An An Gn An Gn An An Tn Tn An An Tn Tn Cn Tn An Tn Tn Tn Gn Gn An An Gn Cn An Tn Tn An Tn An Tn Cn An Gn Tn An Tn Cn Cn An Tn Gn An Cn Tn Gn Tn An An An An Tn Tn Cn An Tn Gn Cn Tn An Tn Gn Tn An Tn An Gn Cn An Tn An Gn An An An Cn An Cn Cn Tn Gn Cn Cn Tn Tn Cn Cn An Tn Tn Cn Tn Tn An Gn An Tn An An Gn Tn Tn Cn Tn Gn Gn An Gn Gn Gn Gn Gn An An An An Gn Gn Cn An Cn Tn An Tn Tn Gn Tn Cn Tn Cn Tn Tn Tn Tn Cn Cn Gn Cn Cn Tn Tn An Gn Tn Cn Tn Cn Cn Tn An An Gn An An Cn Cn Cn Cn An Tn Tn Tn Tn Tn An Cn Tn Tn An Gn Cn Cn An An An Cn Cn An An Cn Tn Gn Cn Cn Tn Cn Tn Gn Gn An Gn An An An An Gn An An An An Gn Tn An Gn Tn Tn Gn An An Gn Cn An Tn Tn Tn Tn Tn Gn Gn Tn Tn An Tn Gn Tn Gn Tn An Gn An Gn Gn Tn Gn Tn An Tn Gn Tn Cn Gn An An An Gn Gn An An An An Gn Cn Cn Cn An Tn Tn Gn An Cn Tn Tn Gn Gn An Gn Cn Tn Tn Gn Tn Tn Tn Tn Gn Cn Tn An Gn Tn Gn An Tn An Tn Cn An Tn Tn Gn An Tn Cn Tn Tn An Cn An Tn An Tn Gn Tn Gn An Gn Tn An Tn Tn Cn An Tn An Tn Cn An An Tn Gn Tn An Gn Tn An Gn Gn An Tn Gn Tn Tn Tn An Tn An Tn Cn Tn Gn Gn An Tn Gn Tn Tn Tn Gn Cn Cn Cn Tn Tn Cn Tn Gn Tn Gn Tn Gn Tn Tn Tn Gn Tn Cn Tn Tn An Gn Gn An Tn Gn Cn Gn Tn Gn Gn Tn Cn Cn Tn An Tn Tn An An Cn An An An Tn Tn Tn Gn Tn Cn An An Gn An Gn Tn Tn Gn Gn An Gn Tn Gn Tn An Tn Tn Tn Gn Tn Tn Tn Tn Cn Tn An Tn Tn Gn Cn Tn Gn Tn Gn Tn An An Cn An An An Tn Tn Tn An Gn Cn An An Cn An An Cn An Tn An Cn An Tn Tn Tn An An Gn An Tn Cn Tn Tn Cn Cn An Gn Tn Tn Tn Cn Cn An Tn Gn Gn Gn Tn Cn An Gn Gn An Cn Tn Cn Cn Cn An Tn Cn An Cn An Gn Cn Tn Tn An Gn Cn Tn Gn An Gn Tn Gn Cn Tn Cn Tn Gn Cn Tn Gn An Gn An Gn Tn Cn Tn Tn Cn Cn An An Gn Gn Cn Tn Gn Cn An Gn Tn Tn An Gn Gn Gn Tn Gn Tn Tn Gn An Cn Tn Tn Tn Gn Cn Tn Gn Tn Gn Tn Tn Tn Cn Tn Tn Tn Tn Tn Gn Cn An Gn Cn Tn Cn Tn Tn Cn Cn An Gn Gn Cn Tn Cn An Cn An An Gn Gn Tn Cn An Cn Tn Gn Gn Cn An Gn An An Tn Tn Cn An Gn Tn Tn Cn Tn Tn Tn Gn Tn Gn Gn Tn Tn Gn Tn An Gn Gn An Cn Tn Gn An Gn Gn Tn Tn Cn Cn Tn Gn Cn Tn Tn Tn Cn Tn Tn Gn Cn Tn Gn Cn Cn Cn An Tn Cn An Gn Tn Gn Gn Gn Gn Gn Gn Cn Cn An Cn Tn Cn Tn Cn An Gn Cn Tn Cn Cn Tn An Gn An Gn Gn Cn Cn Cn Tn Tn Gn Cn Cn Tn Tn Gn Tn Gn Gn Cn An Cn Tn Cn Tn Tn Gn Tn An Gn Gn Cn Cn Cn Tn Cn Tn Cn Cn Tn Gn Gn Cn An Tn Gn Gn Cn An Gn Cn Tn Tn An Cn Cn Tn Cn Tn Tn Cn An An Gn Tn Tn Cn An Gn Gn An Gn Gn An Gn An Cn Tn Tn Tn Cn Cn An Tn Cn Cn An Gn An Cn Tn Gn Cn Tn An An Gn An Tn An Gn An Gn Tn Cn Tn Tn An Cn An Tn An An Tn Gn Tn An Gn Cn An An Tn Cn Cn Cn An An Gn An Gn Tn Gn An Cn An Tn Tn Tn Cn An Tn Cn An Cn Cn Cn Tn Tn Gn Tn Cn An Tn An Tn Tn Cn Tn Gn Tn Tn Gn Gn Cn Tn An Gn An An Gn Cn An An Gn Tn Cn An Tn An Gn Gn Tn Tn Cn Cn An Cn Cn Cn An Cn An Cn Tn Cn An An Gn Gn Gn Gn An Tn Tn An Tn Gn Cn An An Gn Gn Gn An Gn Tn Gn An Gn Tn Cn An Tn Tn Gn Gn Gn An Tn Cn An Tn Tn Tn Tn An Gn An An Tn Tn Cn Tn Gn Cn Cn Tn An Cn Cn An Tn Gn Cn Tn Gn Tn An Tn Gn An An Cn Tn Gn Tn Tn Tn Tn Tn Gn Tn Tn Tn Cn An An Cn Cn An An Cn Tn Tn Cn Tn Gn An Gn An Tn Cn Tn Gn Gn An An An Gn Gn Gn An An Gn An Tn An An Gn An An An Gn Gn Gn Cn Tn Gn Tn Gn Tn An Cn Gn Cn An Gn Tn Tn Gn Cn Cn An Gn Tn An An An Gn An Tn Tn Gn Cn Tn Tn Tn Tn Cn Tn Tn Tn Cn Cn An Cn Cn An An An Gn Cn Tn An An An Gn An Tn An Tn Tn Tn Gn An Gn Gn An An Cn Tn Gn Gn An Tn An Tn An Gn Gn Cn Cn Tn Gn An Gn An Tn Tn Tn Gn Gn Cn An Gn An Gn Cn Tn Cn Tn Gn An Gn Gn Cn Tn Gn Gn An Tn Gn Tn Tn Gn Cn An An Gn Gn An Gn An An An An An An Tn An Gn Gn An An Gn Cn Cn Cn An Cn An Gn Gn Gn Cn Cn An An Gn Cn Tn Tn Gn Gn Gn Cn Cn Tn Cn Cn Tn Tn Gn Tn An Cn Cn Tn Cn Cn Tn Cn Cn An Cn Tn Gn Tn An Tn An Cn An Tn An Cn An Tn Tn Tn Tn Gn Gn Gn Tn Tn Cn An Tn An Tn Tn Tn Tn Tn Cn An Gn Ga Ca Ta Ga Ga Ca Ta Aa Ca Ta Aa Ca Ta Ga Ca Aa Ga Ga Ta Ca Ta Ca Ca Aa Ga Ca Aa Ga Aa Ga Ta Ca Ta Ta Ca Aa Ca Ga Ga Aa Aa Ga Ga Aa Ga Aa Aa Ca Ca Ta Ca Ta Ga Ga Ca Ca Ta Ta Ga Aa Ga Ga Ta Ga Ta Ca Aa Ta Ga Ca Ga Ta Ga Ga Aa Aa Ga Ga Aa Ta Aa Aa Ga Ca Ta Ga Ga Ta Ga Ta Aa Ca Aa Aa Ta Ga Ta Ga Ca Ta Ta Ta Aa Ca Ta Aa Ta Ca Ga Aa Aa Aa Ta Ga Ga Ca Aa Aa Aa Ga Ca Ca Ta Ta Ta Aa Aa Ga Ta Ta Ta Ta Ta Ca Ca Aa Ca Ta Ga Ga Aa Aa Ta Ta Ca Ta Aa Aa Ca Ca Ta Ca Aa Ca Ca Aa Ta Ta Ca Ta Ga Aa Aa Aa Aa Ca Ca Aa Aa Ca Aa Ta Aa Aa Ga Ta Ca Aa Ca Aa Aa Ta Ga Ga Ca Aa Ca Ca Ta Aa Ca Ca Aa Ta Ta Ga Ca Ta Ca Aa Ga Ga Ca Aa Ta Ga Ga Ga Aa Aa Aa Ga Ca Aa Ta Ca Ga Ca Ta Aa Ca Aa Ca Aa Ta Ca Aa Ga Ca Aa Ga Ga Aa Aa Ta Aa Ta Ca Ta Ga Ta Ca Aa Ca Ta Ga Ta Ga Aa Aa Aa Ga Gn Tn An Tn Tn Gn Tn An Tn Tn Gn Gn An An Tn An Gn Tn Cn An Tn An Gn An An Cn Tn Gn An Tn An Gn Tn Cn Cn Cn Tn Cn Cn Cn Cn Cn Tn Gn An Gn Gn Gn An Cn Cn An Tn Cn An Tn An An An Tn An Tn Tn Cn Tn An An An Cn Tn Cn Cn Tn Cn An Cn Tn Tn Tn An An Tn Tn Tn An Cn An An Gn Tn Gn An An Gn An An An Cn Cn Gn Gn Gn Cn Tn Cn Cn An Gn An Gn An Gn Gn Tn An An An Cn Tn Gn Cn An Tn Tn An Cn Tn An An An Gn Gn Cn Cn An Cn An Cn Cn An Tn Gn An Cn Cn An Gn Tn An Gn Cn Tn Gn Gn An An Cn Cn An Gn An An Cn Tn Gn An Gn Gn Tn Cn Tn Cn Tn Gn Gn Gn Tn Cn Cn Tn An Gn Tn Cn Cn An An Tn An Cn Tn Cn Tn Tn Tn Cn Cn An Gn Tn Gn Tn An Cn Cn Gn Cn An An Cn Cn Cn Cn Cn An Tn Cn An An Cn An An Tn Cn An Cn An Gn Gn An Gn Cn Tn An An Tn Gn Tn Cn An An Tn Gn Tn Cn An Gn Gn Gn Cn An Tn An An Tn Gn Gn Tn Gn An Tn Cn Cn An Tn Tn Tn Cn An Tn Cn Cn An Tn Cn An Tn Tn An Tn Tn Tn Tn An Gn Tn Tn Gn Tn An An An Gn An An Tn An Cn An Gn An Cn Tn Cn An Cn Tn Cn Tn An Gn Gn Tn An Gn An Tn An Gn Gn An Gn An Gn Tn Tn An Gn An An An An Tn An An Tn An Cn Cn An Tn Gn Gn An An Tn Cn Tn Cn An Tn Gn Gn An An Gn Cn Cn Cn Cn An Tn An Gn An An An Gn An Gn An An An Cn Cn Tn Gn Gn An An An An Tn An Tn Gn An An Gn An An Cn Gn Gn Gn Cn Tn An Tn An Cn Tn Gn Tn Cn Cn An Tn Cn Tn An Cn Cn Tn Cn Tn Cn An Gn Gn Gn An Cn An An An Gn Cn An Gn Cn Cn Tn Cn Tn Gn Tn Gn Gn Tn Gn Cn Cn An Cn Tn Tn Cn Tn Tn Cn Cn Tn Gn An Gn Cn An Tn Cn Tn Gn Cn An Gn Gn Cn Tn An An Gn Cn An Gn Tn Tn Tn Gn Cn An Tn An Tn Gn Gn Cn Cn Tn Tn Gn Tn Tn Cn Tn Gn An Cn Cn Cn Tn Tn Tn Cn An Cn Gn Gn Cn Cn Tn Tn Tn Cn Cn Tn Cn Tn Tn An An Tn An Cn An An Tn Cn An Cn Tn An Gn Cn An An Cn Cn Tn Tn Tn Cn Tn Cn Tn Gn Tn Cn Tn Gn An Gn An An An An An An An An An An Gn Tn Tn Cn An Cn Tn Cn Cn Tn Gn Tn Gn Gn Gn An An An Gn An An An Tn Tn An Tn Tn Gn Tn Cn Cn Cn An Gn An Tn Tn An Tn Cn Tn Tn Tn Tn An An An Gn Tn Cn An Gn Gn An An Cn An An An An An Tn An An Tn Gn An Gn Tn Cn An Cn Tn Gn An An An Tn Gn Tn Gn Gn Cn Gn Gn Gn Cn An Cn Cn Tn Gn Tn An Gn Tn Cn Cn Cn An Gn Cn Tn An Cn Tn Cn Gn Gn Gn An Gn Gn Cn Tn Gn An Gn Gn Cn An Gn Gn An Gn An An Tn An Gn Cn An Tn Gn An An Cn Cn Cn Gn Gn Gn An An Gn Cn Gn Gn An Gn Cn Tn Tn Gn Cn An Gn Tn Gn An Gn Cn Tn Gn An Gn An Tn An Gn Cn Gn Cn Cn An Cn Tn Gn Cn An Cn Tn Cn Cn An Gn Cn Cn Tn Gn Gn Gn Cn An An Cn An Gn An Gn Cn Gn An Gn An Cn Tn An Cn Gn Tn Cn Tn Cn An An An Tn An An Tn An An Tn An An Tn An An Tn An An Tn An An Tn An An Tn An An Tn Gn An Tn Gn An Tn An An Tn An An An Gn Tn Cn An Cn Tn Gn An An An Tn Gn Gn Cn An An Gn Tn Gn Tn Gn Tn Gn An Gn Tn Cn An Gn Gn Gn Tn Tn Cn An Cn Tn Cn Cn Tn An Gn Tn Tn Gn An Cn Cn An Cn Tn Gn Cn Cn An An An An Gn Gn Gn Cn An Tn Gn Gn Tn Tn An Cn Tn Gn An Cn Cn An Gn An Tn Gn An Cn An Cn Gn Gn Tn Cn An Tn Tn Cn An Gn Gn Gn An An Gn Gn An Tn Gn Gn An Cn Cn Tn Tn Gn Tn An Gn An Gn Cn An Gn Gn Gn Gn Tn Tn Cn Cn Cn An An Cn Cn Cn Tn Tn An Cn Tn Gn Gn Tn Cn Cn Gn Tn Gn Gn Cn Cn Tn Gn Tn Tn An Gn Gn An An Cn Tn Gn Gn Gn Cn Tn Gn Cn An Cn An Gn Cn An Gn Gn An Gn Gn Tn Gn An Gn Cn Tn Tn Cn Tn Tn Cn An Tn Gn An Gn Cn Cn An Gn Cn An Tn Tn An Cn Gn Gn Cn Cn Tn Gn An Gn Cn Tn Cn Cn Tn Cn Cn Tn Cn Cn Tn Gn Tn Cn Cn An An Tn Cn An Gn Tn Gn Gn Cn An An Cn An Tn Tn An Gn An Tn Tn Cn Tn Cn An Cn An Gn An An An Cn An Tn Gn An An Cn Cn Cn Tn An Tn Tn Gn Tn Gn An An Tn Tn Gn Tn Gn Cn Gn Tn Gn Cn Cn Tn Gn An Tn Gn An Tn Cn Tn Gn An Gn Gn Tn Gn Gn An An Cn An Gn Tn Tn An Cn Cn Tn Tn Cn Cn An An An An Cn Tn Gn Tn Cn Cn Cn Cn An Cn Tn Tn Cn An Cn Cn Cn Cn Cn Cn Gn Gn Cn Tn Gn Tn Gn Gn An An An An Gn Tn Tn Gn Cn Cn Tn Tn Cn Cn An Cn An An An An Tn Cn Cn An Tn Cn Cn Cn Tn An Gn Tn Gn Cn Cn An An An An An Gn Gn Tn Tn Gn Gn Gn Gn An Cn Cn An Cn Tn Gn Gn Tn An Tn An Gn An Gn Gn Tn An Tn Cn Cn Tn Cn An Gn Tn An Gn Gn An Cn Tn Cn An Gn Cn An An Gn Cn Tn Gn Gn Cn An An Cn Cn Cn Tn An An Tn Tn An Tn Gn Tn Cn Tn An Tn Tn An Gn Gn An Cn An Cn Cn Cn Cn An An Gn An An Tn Gn Gn Cn Tn Cn Tn Cn Tn Gn Cn Tn Gn Gn An An Gn Tn An An An An An Gn Gn Gn Tn Tn An An Tn Gn Cn Cn Tn Tn An Tn Gn Gn An Tn Gn Cn An Tn Tn Tn Tn Cn Tn An Gn Gn Gn Cn Cn An Gn Gn Tn Tn Tn Cn Cn Cn An An Gn Gn Gn Gn Gn Tn An An An Gn Gn Gn Cn An Tn Gn Tn Cn Tn Tn Tn Tn Gn Tn Gn An An An An Gn Gn An Cn Cn Tn Gn Gn An Tn Gn Cn Tn An An An Cn An Gn Gn Cn An An Cn Cn Tn Tn Gn Tn Cn Cn Cn Cn Cn An Tn Cn An An Cn Tn Tn Tn Cn Tn Cn Cn Tn Tn An Gn Aa Ga Ca Ta Aa Ta Ta Ta Ca Ca Aa Ga Ca Ta Ca Ca Aa Ga Ta Ga Ca Ta Ga Aa Aa Ta Ga Ca Aa Ta Ca Ta Ga Ta Ga Aa Ca Aa Ta Ca Ca Ca Ca Aa Ca Ta Ca Ca Ta Ga Ga Aa Ga Ga Ga Ga Aa Aa Ta Ca Ta Ga Ga Ta Ca Aa Ca Ca Ca Ta Ga Aa Ga Ca Ta Ga Ta Ga Aa Aa Aa Ca Aa Aa Aa Ga Ta Ta Ga Ca Ta Ca Ta Ta Ga Ca Aa Ga Aa Ga Ga Ca Ca Ta Ga Ga Ta Ta Ta Ga Ca Aa Ga Ca Ta Ta Ta Aa Ca Ta Ta Ca Ta Ca Ca Ta Ta Ca Ta Aa Ca Aa Ta Ga Ga Ga Ca Aa Ga Ca Aa Aa Ga Aa Ca Ca Ca Ta Ga Ca Ga Aa Ga Ga Ca Aa Ga Ga Aa Aa Ca Aa Ca Aa Ta Ca Ca Ta Ca Ta Ga Aa Aa Ta Aa Ca Ca Aa Aa Aa Ta Aa Ca Ta Aa Aa Ca Ta Ga Ca Ta Aa Ga Aa Aa Ga Aa Ga Aa Aa Ga Aa Ca Ta Ca Ta Ga Ga Ga Ta Ta Aa Ta Aa Ca Ta Ga Ga Ta Ga Ca Ga Aa Ga Ga Ca Ta Ga Ca Ca Aa Ca Aa Ga Aa Ga Ga Aa Ta Ga Ga Aa Aa Aa Ta Ga Ta Ca Ca Ta Ta Aa Aa Ga Ca Ga Ca Aa Ga Ca Ca Ca Ta Ga Aa Ga Ta Ta Ga Ga Aa Ga Ca Ta Ta Ca Aa Aa Ga Ta Ga Ca Ta Ta Ga Gn Tn Gn An Gn Tn Gn An Gn An An Tn Gn An Cn Gn Gn Gn An An Gn Cn Cn An Cn Tn Gn Gn Cn An Cn An Gn An An Gn An An Gn Gn Gn An Cn Tn Cn Cn Cn Tn Tn An Tn Cn Tn Cn Cn Cn An Tn Gn Gn Gn An Cn Tn Gn An Gn Gn Tn Tn Tn Gn Tn Tn Cn An An Gn Gn Gn Tn Tn Tn Tn Tn Gn Gn Cn Cn Cn An Gn An Cn An Gn Gn An Gn Gn Gn Gn An An An Gn Tn Cn Tn Cn Tn Tn Cn An Gn Gn An An An An Gn Cn Cn Cn An Cn An An Gn Cn An Gn Gn Cn Cn Tn Tn Tn Cn Cn An Tn Cn Cn Tn Tn Gn An Tn Tn Cn An Cn An An Cn An Tn Cn An Cn Tn Cn Tn Tn Cn Tn Cn Cn Tn Cn Gn Cn An An An Cn Tn Gn Tn Tn An An An Tn Tn Tn Cn Cn Tn Tn Tn Cn Cn Tn Tn Tn Cn Tn Tn Tn Tn Tn Cn Tn Tn Tn Tn Tn Cn Cn Tn Tn Tn Gn Cn Cn Tn Tn Tn Cn Cn Tn Tn Cn Cn Tn Cn Cn An Tn Tn Tn Cn Tn Tn Tn Cn Cn Tn Tn Cn An Tn Tn Tn Tn Cn Tn Cn Cn Tn Cn Tn Gn Tn Cn Cn Tn Tn Cn Tn Tn Tn Tn Tn Tn Tn Cn Tn Cn Cn Tn Tn Cn An Tn Tn Tn Tn An Tn Tn Tn Tn Cn Cn Cn Cn Tn Cn Cn Cn Tn Cn Cn Cn An Cn Tn Cn Tn Tn Cn Cn Cn Tn Cn Cn An Cn Tn Cn Cn An Tn Gn An Cn Cn Cn Cn Cn Cn Cn Tn Tn Cn Tn Cn Tn Cn Tn Cn Tn Cn Tn Cn Tn Cn Tn Cn Cn Cn Tn Gn Cn Tn Tn Cn Cn Cn Tn Gn Cn Cn Tn Cn Cn Cn Tn Cn Cn Tn Cn Cn Tn Cn An Cn Cn An An Cn An An Tn Cn Tn Cn An Cn Cn An An Cn An An Tn Tn Tn An Tn Cn An An Gn Tn Tn Cn Tn Tn Cn Cn Tn An Tn Cn Tn Gn Tn Tn Gn Tn Cn An Tn An Cn Gn Tn Cn Tn Gn Gn Gn Gn An Tn An Tn An An An Gn An Cn An Tn Tn Tn Gn An Gn Tn An Tn An Gn Tn Tn Cn Tn Tn Gn Cn Tn Tn Cn An An Gn Gn An Gn Cn Tn Cn An Cn An An Gn Gn Tn Gn Gn An Tn Tn Tn An Tn Cn An Gn An Cn An Gn Tn Gn An Tn Tn Tn Tn Gn Tn An An An Cn Tn Gn Cn An An An Tn Cn An Cn Cn An Cn Cn Tn Cn Cn Cn Cn An An Gn Tn An Tn Cn Tn Cn Tn An Tn Tn Tn An An Cn Tn Gn An Gn Gn Cn An Gn An Gn Gn Gn An Tn Tn Gn Tn Gn An Gn Cn Tn Cn Tn An Gn An An Cn An An An Gn Tn Cn Tn Cn Cn Tn Tn Gn Gn Gn Gn Gn Gn Gn An An An An An An Gn Tn Tn Cn An Tn Cn Tn Tn Cn An An Cn Cn Cn An An An Tn Tn Cn An Tn Tn Tn Cn An An Gn Tn An Tn Tn An An An Tn Gn Gn Cn An Cn An Gn An Gn An Tn An Tn Cn An Gn Tn Tn Gn Tn Cn Tn Cn Tn Gn Gn An An Cn Tn An Gn Gn Gn An Gn Tn An An Gn Tn Cn Cn An Cn Tn Gn An Cn An Gn Gn Gn Cn Cn Cn An Gn Cn An An Tn Tn An An Gn Cn Tn Cn Tn Tn Cn Cn An An Gn Gn An Gn Cn Cn Tn Gn Tn Cn Cn Cn Tn Gn Tn Cn Tn An Tn An Cn An Cn Cn An Tn An Tn An Gn Cn Cn An Gn Tn Tn An Gn An Gn Cn Tn An Cn Cn Gn Cn Cn An Gn Tn Tn Cn Cn Tn Tn Cn Cn Tn Gn Cn Tn Cn Cn Cn Tn Gn An An An Tn Gn An Cn Cn An Gn Tn Gn Cn Cn Tn Cn Cn Cn Tn Gn An Gn Gn An Cn An Gn Gn Cn An Cn An Tn Cn An Gn Gn Tn Cn Tn Cn An Gn Cn Cn An An Cn Cn Tn Tn Cn Cn Cn Tn Tn Cn Cn An An An Cn An Tn An An Cn Tn Cn An Gn Cn Tn An Gn An Cn Cn Cn Cn Tn Cn Tn Gn Gn Tn Cn Tn Cn Tn An An An Tn An Gn Tn An Tn Cn Tn Cn Tn Tn Cn Tn Cn Tn Tn Tn Gn Tn Cn Tn Tn Tn Tn Tn Cn Tn Gn Tn Tn Tn Cn An Gn Ga Ca Ca Ta Ca Ca Aa Ga Ta Ta Aa Ca Ca Aa Aa Ca Ta Ca Ca Ta Ga Ta Ca Ta Ga Ga Ta Ta Ta Ca Aa Ta Ga Ta Ca Ca Ta Ta Ta Ta Ca Ta Aa Ta Ca Ta Ga Ga Ca Aa Ga Ta Ga Ga Ga Aa Aa Ta Aa Aa Ta Ga Ta Ta Ta Ta Ta Aa Ga Ta Ga Aa Aa Ca Aa Ca Ta Ga Ta Ta Ca Ta Ca Ta Ga Ga Ga Ta Ga Aa Ca Aa Aa Ta Aa Ca Ga Ta Aa Aa Aa Ga Aa Aa Ca Ta Ga Aa Aa Aa Aa Ga Aa Aa Aa Ga Aa Aa Aa Aa Aa Ga Ta Ga Ga Ga Aa Ta Ta Ta Aa Ga Aa Aa Aa Ta Ca Ta Ca Ta Ta Ta Ga Ga Aa Ta Ta Ca Ta Ga Ga Ta Ca Aa Tf Gf Af Gf Af Af Gf Af Af Gf Gf Tf Af Af Tf Tf Tf Cf Cf Af Gf Cf Cf Tf Tf"
    #seq, a, string = get_testable_string(anno)
    #print(seq, len(seq))
    #print(a, len(a))
    #print(string)

    #path = predict_all(seq, string)
    #calculate_accuracy(a, path)
    #lines = get_annotated_lines('CDS', route, 'complement')

    #with open('out_CDS_complement_annotated_gene.txt', 'w') as o_file:
    #    for l in lines:
    #        for line in l:
    #            if line:
    #                o_file.write(line + '\n')

    with open('out_CDS_complement_annotated_gene.txt') as i_file:
        uniques = set()
        counts = 0

        sensitivity_accumulator = 0
        specificity_accumulator = 0
        cc_accumulator = 0
        average_specificity = 0
        average_sensitivity = 0
        true_positive_a = 0
        true_negative_a = 0
        fake_positive_a = 0
        fake_negative_a = 0

        for i, line in enumerate(i_file):
            parts = line.split('-')
            gene_id = parts[1].replace('\n', '')
            anno = parts[0]
            if gene_id not in uniques:
                uniques.add(gene_id)
                seq, a, string = get_testable_string(anno)
                path = predict_all(seq, string)
                res = calculate_accuracy(a, path)
                if res:
                    print(gene_id, res)
                    sensitivity_accumulator += res['sensitivity']
                    specificity_accumulator += res['specificity']
                    cc_accumulator += res['cc']
                    true_positive_a += res['true_positives']
                    true_negative_a += res['true_negatives']
                    fake_positive_a += res['fake_positives']
                    fake_negative_a += res['fake_negatives']
                    counts += 1

        average_sensitivity = sensitivity_accumulator / counts
        average_specificity = specificity_accumulator / counts
        average_cc = cc_accumulator / counts
        print(average_sensitivity, average_specificity, average_cc, true_positive_a, true_negative_a, fake_positive_a, fake_negative_a)
