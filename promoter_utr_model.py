from matrix_from_aln import matrix_from_fasta
import numpy
from model_maker_utils import sequence_state_factory
from model_maker_utils import classify
from model_maker_utils import add_sequence
from model_maker_utils import spacer_states_maker
from model_maker_utils import add_variable_length_sequence
from model_maker_utils import load_long_training_examples
from converter_to import converter_to
import calculator
from pomegranate import State
from pomegranate import HiddenMarkovModel
from pomegranate import DiscreteDistribution


matrix_TATA = numpy.array(matrix_from_fasta('tata_-5_11_completa.seq'))
matrix_GC = numpy.array(matrix_from_fasta('gc_completo.seq'))
matrix_CCAAT = numpy.array(matrix_from_fasta('CCAAT_completa.seq'))
matrix_Inr = numpy.array(matrix_from_fasta('Inr_completo.seq'))
matrix_no_inr = numpy.array(matrix_from_fasta('no_inr.fa'))

gc_data = classify(matrix_GC, 2)
tata_data = classify(matrix_TATA, 2)
cat_data = classify(matrix_CCAAT, 2)
inr_data = classify(matrix_Inr, 2)
no_inr_data = classify(matrix_no_inr, 2)

no_coding = calculator.intron_calculator('cuts_intron.txt')


# Model
promoter_utr_model = HiddenMarkovModel('promoter')

# States
back = State(DiscreteDistribution(no_coding.p), name='back')

gc_states = sequence_state_factory(gc_data, 'GC')
post_gc_var_spacers_tss = spacer_states_maker(151, no_coding.p, 'post gc var spacer tss')
post_gc_spacers_tss = spacer_states_maker(38, no_coding.p, 'post gc spacer tss')

post_gc_var_spacers_tata = spacer_states_maker(151, no_coding.p, 'post gc var spacer tata')
post_gc_spacers_tata = spacer_states_maker(18, no_coding.p, 'post gc spacer tata')


cat_states = sequence_state_factory(cat_data, 'CAT')
post_cat_var_spacers_tss = spacer_states_maker(151, no_coding.p, 'post cat var spacer tss')
post_cat_spacers_tss = spacer_states_maker(42, no_coding.p, 'post cat spacer tss')

post_cat_var_spacers_tata = spacer_states_maker(151, no_coding.p, 'post cat var spacer tata')
post_cat_spacers_tata = spacer_states_maker(22, no_coding.p, 'post cat spacer tata')

tata_states = sequence_state_factory(tata_data, 'tata')
post_tata_var_spacers = spacer_states_maker(16, no_coding.p, 'post_tata_var_spacer')
post_tata_spacers = spacer_states_maker(4, no_coding.p, 'post_tata_spacer')

inr_states = sequence_state_factory(inr_data, 'inr')

no_inr_states = sequence_state_factory(no_inr_data, 'no inr')
# Add States
promoter_utr_model.add_state(back)


# Add Sequences

#GC

add_sequence(promoter_utr_model, gc_states)

add_sequence(promoter_utr_model, post_gc_spacers_tata)
add_variable_length_sequence(promoter_utr_model, post_gc_var_spacers_tata, post_gc_spacers_tata[0])

add_sequence(promoter_utr_model, post_gc_spacers_tss)
add_variable_length_sequence(promoter_utr_model, post_gc_var_spacers_tss, post_gc_spacers_tss[0])

add_sequence(promoter_utr_model, inr_states)
add_sequence(promoter_utr_model, no_inr_states)

# CAAT
add_sequence(promoter_utr_model, cat_states)

add_sequence(promoter_utr_model, post_cat_spacers_tss)
add_variable_length_sequence(promoter_utr_model, post_cat_var_spacers_tss, post_cat_spacers_tss[0])

add_sequence(promoter_utr_model, post_cat_spacers_tata)
add_variable_length_sequence(promoter_utr_model, post_cat_var_spacers_tata, post_cat_spacers_tata[0])


# TATA
add_sequence(promoter_utr_model, tata_states)
add_sequence(promoter_utr_model, post_tata_spacers)
add_variable_length_sequence(promoter_utr_model, post_tata_var_spacers, post_tata_spacers[0])


# Transitions
promoter_utr_model.add_transition(promoter_utr_model.start, back, 1)

promoter_utr_model.add_transition(back, back, 0.99)
promoter_utr_model.add_transition(back, tata_states[0], 0.00043)
promoter_utr_model.add_transition(back, gc_states[0], 0.00707)
promoter_utr_model.add_transition(back, cat_states[0], 0.0024)
promoter_utr_model.add_transition(back, inr_states[0],  0.000031)
promoter_utr_model.add_transition(back, no_inr_states[0],  0.000069)

promoter_utr_model.add_transition(gc_states[-1], post_gc_var_spacers_tata[0], 0.1)
promoter_utr_model.add_transition(gc_states[-1], post_gc_var_spacers_tss[0], 0.9)

promoter_utr_model.add_transition(post_gc_spacers_tata[-1], tata_states[0], 1)

promoter_utr_model.add_transition(post_gc_spacers_tss[-1], inr_states[0], 0.35)
promoter_utr_model.add_transition(post_gc_spacers_tss[-1], no_inr_states[0], 0.65)

promoter_utr_model.add_transition(cat_states[-1], post_cat_var_spacers_tss[0], 0.86)
promoter_utr_model.add_transition(cat_states[-1], post_cat_var_spacers_tata[0], 0.14)

promoter_utr_model.add_transition(post_cat_spacers_tata[-1], tata_states[0], 1)

promoter_utr_model.add_transition(post_cat_spacers_tss[-1], inr_states[0], 0.32)
promoter_utr_model.add_transition(post_cat_spacers_tss[-1], no_inr_states[0], 0.68)


# TATA
promoter_utr_model.add_transition(tata_states[-1], post_tata_var_spacers[0], 1)

promoter_utr_model.add_transition(post_tata_spacers[-1], inr_states[0], 0.31)
promoter_utr_model.add_transition(post_tata_spacers[-1], no_inr_states[0], 0.69)

promoter_utr_model.add_transition(inr_states[-1], promoter_utr_model.end, 1)
promoter_utr_model.add_transition(no_inr_states[-1], promoter_utr_model.end, 1)


promoter_utr_model.bake()

with open('promoter_utr_model_base.json', 'w',  encoding='utf-8') as out:
    out.write(promoter_utr_model.to_json())

string = """CATTCCCAGTTCTTTACATTCATCCCTTGTTTCCAGAAAGGGCAGAGGAAGCGAGGAAAA
AGTGCGTGGCCTGAAGTGACGCCTGGCGTTGCCCGAAGCCCGCCCAGCGCTGCCAGGTGA
CGCCACTGCGACACAAAGGCGGGGATTGCGTAGGGAAAGGCCCTAGGCCATAAACGGGGG
TGGGGCCTCCCCGGAGGCCAGTGCG"""
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

# promoter_utr_model.fit(load_long_training_examples('representative_promo_-499_5.fa', 1002), n_jobs=4, max_iterations=100)

# with open('promoter_utr_model_trained_representative.json', 'w',  encoding='utf-8') as out:
#     out.write(promoter_utr_model.to_json())
