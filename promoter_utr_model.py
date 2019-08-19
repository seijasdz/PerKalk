from matrix_from_aln import matrix_from_fasta
import numpy
from model_maker_utils import sequence_state_factory
from model_maker_utils import classify
from model_maker_utils import add_sequence
from model_maker_utils import spacer_states_maker
from model_maker_utils import add_variable_length_sequence
from converter_to import converter_to
import calculator
from pomegranate import State
from pomegranate import HiddenMarkovModel
from pomegranate import DiscreteDistribution


matrix_TATA = numpy.array(matrix_from_fasta('tata_-5_11_completa.seq'))
matrix_GC = numpy.array(matrix_from_fasta('gc_completo.seq'))
matrix_CCAAT = numpy.array(matrix_from_fasta('CCAAT_completa.seq'))

gc_data = classify(matrix_GC, 2)
tata_data = classify(matrix_TATA, 2)
cat_data =classify(matrix_CCAAT, 2)

no_coding = calculator.intron_calculator('cuts_intron.txt')


# Model
promoter_utr_model = HiddenMarkovModel()

# States
back = State(DiscreteDistribution(no_coding.p), name='back')

gc_states = sequence_state_factory(gc_data, 'GC')
post_gc = State(DiscreteDistribution(no_coding.p), name='post_gc')
post_gc_spacers_tss = spacer_states_maker(150, no_coding.p, 'post_gc_spacer_tss')
post_gc_spacers_tata = spacer_states_maker(70, no_coding.p, 'post_gc_spacer_tata')

cat_states = sequence_state_factory(cat_data, 'CAT')
post_cat = State(DiscreteDistribution(no_coding.p), name='post_cat')
post_cat_spacers_tss = spacer_states_maker(80, no_coding.p, 'post_cat_spacer')
post_cat_spacers_tata = spacer_states_maker(30, no_coding.p, 'post cat spacer tata')

tata_states = sequence_state_factory(tata_data, 'tata')
post_tata_var_spacers = spacer_states_maker(16, no_coding.p, 'post_tata_var_spacer')
post_tata_spacers = spacer_states_maker(7, no_coding.p, 'post_tata_spacer')

# Add States
promoter_utr_model.add_state(back)
promoter_utr_model.add_state(post_gc)
promoter_utr_model.add_state(post_cat)

# Add Sequences
add_sequence(promoter_utr_model, gc_states)
add_sequence(promoter_utr_model, tata_states)
add_sequence(promoter_utr_model, cat_states)
add_sequence(promoter_utr_model, post_tata_spacers)
add_sequence(promoter_utr_model, post_gc_spacers_tss)
add_sequence(promoter_utr_model, post_gc_spacers_tata)
add_sequence(promoter_utr_model, post_cat_spacers_tss)
add_sequence(promoter_utr_model, post_cat_spacers_tata)

add_variable_length_sequence(promoter_utr_model, post_tata_var_spacers, post_tata_spacers[0])
# Transitions
promoter_utr_model.add_transition(promoter_utr_model.start, back, 1)

promoter_utr_model.add_transition(back, back, 0.994)
promoter_utr_model.add_transition(back, tata_states[0], 0.001)
promoter_utr_model.add_transition(back, gc_states[0], 0.003)
promoter_utr_model.add_transition(back, cat_states[0], 0.002)

promoter_utr_model.add_transition(gc_states[-1], post_gc, 1)

promoter_utr_model.add_transition(post_gc, post_gc, 0.5)
promoter_utr_model.add_transition(post_gc, post_gc_spacers_tata[0], 0.17)
promoter_utr_model.add_transition(post_gc, post_gc_spacers_tss[0], 0.18)
promoter_utr_model.add_transition(post_gc, cat_states[0], 0.15)  # TODO cambiar por un spacer correspondiente

promoter_utr_model.add_transition(cat_states[-1], post_cat, 1)

promoter_utr_model.add_transition(post_cat, post_cat, 0.5)
promoter_utr_model.add_transition(post_cat, post_cat_spacers_tss[0], 0.431)
promoter_utr_model.add_transition(post_cat, post_cat_spacers_tata[0], 0.069)  # TODO cambiar por un spacer correspondiente
promoter_utr_model.add_transition(post_cat_spacers_tata[-1], tata_states[0], 1)

promoter_utr_model.add_transition(post_cat_spacers_tss[-1], promoter_utr_model.end, 1)

promoter_utr_model.add_transition(post_gc_spacers_tata[-1], tata_states[0], 1)

promoter_utr_model.add_transition(post_gc_spacers_tss[-1], promoter_utr_model.end, 1)

promoter_utr_model.add_transition(tata_states[-1], post_tata_var_spacers[0], 1)

promoter_utr_model.add_transition(post_tata_spacers[-1], promoter_utr_model.end, 1)


promoter_utr_model.bake()

print(promoter_utr_model.to_json())

string = """GCAGTCACTTTTCTTCTACAGCCAGAGGCTATGATGTCAGGAAGTAAATGATAGAAGATG
TCAGCTCTCCTATCTATAACCTCAGCTTTGCAGTTTCTCTGGGGGCCTCATGCACTCTTA
CAGCTTCTTAAACAGGAAGTAGACCTGGCAAAAAATTTAAGTGCGTTTTACACTTGATAT
ACCACAGACTCGGCCAGAACAGTTCTAAACAGAAGACTTCACAGAAATCCTGACACTTTG
TTTTTTGCCAAAGAACAGTGTCTTTTTCAGAAATGCTGATTTTGGGATCATGGACATTCT
CCTTAATTTCTCCTTAGCAAGTCTGTTGGAAGGAGGAAATTGGTGCCTTTGAGCAGATTA
TCCCCCCCACACAACAGCGCTTGATAAGGGACGGTAAGAAGCACAGGTGAGTGCTGTTAG
AAGGTGGGGGTGGGAAAGAGGAGGAAGGCACAGGAGGAGTCAATAAGAGGTTAAATAGAA
AGACCAGAGCAACCCGAGA"""
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


