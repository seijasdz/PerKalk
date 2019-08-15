from matrix_from_aln import matrix_from_fasta
import numpy
from model_maker_utils import sequence_state_factory
from model_maker_utils import classify
from model_maker_utils import add_sequence
from converter_to import converter_to
import calculator
from pomegranate import State
from pomegranate import HiddenMarkovModel
from pomegranate import DiscreteDistribution


matrix_TATA = numpy.array(matrix_from_fasta('tatas-31-18.fa'))

tata_data = classify(matrix_TATA, 2)
intron_distribution = calculator.intron_calculator('cuts_intron.txt')

# Model
tata_model = HiddenMarkovModel()

# States
back = State(DiscreteDistribution(intron_distribution.p), name='back')
post_tata = State(DiscreteDistribution(intron_distribution.p), name='post')
tata_states = sequence_state_factory(tata_data, 'zona tata')

print(tata_states)

# Add States
tata_model.add_state(back)
tata_model.add_state(post_tata)

# Add Sequences States and transitions
add_sequence(tata_model, tata_states)

# Transitions
tata_model.add_transition(tata_model.start, back, 1)
tata_model.add_transition(back, back, 0.9999)
tata_model.add_transition(back, tata_states[0], 0.0001)
tata_model.add_transition(tata_states[-1], post_tata, 1)
tata_model.add_transition(post_tata, post_tata, 1)

tata_model.bake()

print(tata_model.to_json())

string = 'GAGGGCTCTGACTCGCCCAAGGCCACACAGCCTTGCGTGGGCCTTTCACATCCCACACAACAGAGGGGCATCCTCAGCCTGGTTGGCAGAGGGCAGGCAGGATAGATGGCAGAGTCTTCTCCGAGGAGAGGGGTTTTGCTCATGGAACCTCCTCCTCCACACTCACAGCCCTGGGCGAGACCTGTGGAGCAGCCGCCAACAGAGTGAGGGAGGGGGCTCGGGGCAGCTGGGGGTGACTTGAGGAAGTCCAGCTGGACTGCGAGGGGCCCCTGGGGACTGCCAGGGAGCCTCAGGACTCCCAGAGGTGCTCCAGGCACAGAGGGAGGAATGGGCCTTCCATCTCCCTCCCCCCTTCACTGCAGAGGCTGGGTCGGGCCAGGCGCCCGGGGAGGAGGCGGTGTCCCTGGCTCCCAGCCCGCCGGTGCAGCGGGGCAGGGCTGGACCAGAAGGGGTGGGGCACCGTGCCTGGTATAAGAGGCAGCCAGGGCACCGAGGCAATGAGCTATCTGCTCAGCTTAATAGCAGGACGCTGGCAACAGGCACTCCCTGCTCCAGTCCAGCCTGCGCGCTCCACCGCCGCTATGGTCTCCGTGCCTACCA'.lower()
lists = list(string)
two = converter_to(lists, 2)
print(two)
seq = numpy.array(two, numpy.unicode_)

logp, path = tata_model.viterbi(seq)

print([p[1].name for p in path])



