from pomegranate import HiddenMarkovModel
from converter_to import converter_to
import numpy

with open('promoter_utr_model_trained.json') as base_model_file:
    model_json = base_model_file.read()

promoter_utr_model = HiddenMarkovModel.from_json(model_json)

string = """TACATAAGAGAAACGTGGCCAACAGACCGAGGTGGGGACGGGGACAGGGACCGGCAATGCA
GGAAATCCGAGTGTCACATCCTCTGCCTCTCATTTGCACACTGCTCCCTCGCTATGCTCA
CCGCTCCCGCCGATCCAGGGACGTGATCCAGGGACTCTGGGAAATGCAAAGCTACACACA
GTGGAGCGGGGGCTGGGGGTGTGTAGACCGCCGGGATTCCGAGTTTCCCGGCACGCCTAG
GAGAGGGAGAGGCAGGCAATGTCAGGGAAATTGGGCAGGCAAGACGCCAGGGACGCCACG
TACTGCCAGGTTCTCAACGAGGTGGAGCCAAAGGGGCAGGCCCCGCGGTGCGCCCGGCGC
TGGGCTCACGGGTTGCTGCACCCGGCCCAGGATCGCGGGCGGTGCAGACTCAGCAGGGGC
GGGTGCAAGGACGAGGCGGGGCCTCTGCGCCCGGCCCTCTTCCCGGACTATAAAGAGAGC
CGCCGGCTTCTGGGCTCCACCACGC"""
string = string.lower().replace('\n', '')
print(len(string))
lists = list(string)
two = converter_to(lists, 2)
print(two)
seq = numpy.array(two, numpy.unicode_)

logp, path = promoter_utr_model.viterbi(seq)

path_names = [p[1].name for p in path]
print(path_names)
count = 0
print([(string[i + 1], name, i - len(path_names) + 1) for i, name in enumerate(path_names) if i + 1 < len(string)])