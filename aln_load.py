import numpy
from pomegranate import State
from pomegranate import DiscreteDistribution
from insertor import HMMWrapper


def matrixer(filename):
    to_matrix = []
    with open(filename) as file_handle:
        for index, line in enumerate(file_handle):
            tokenize = line.split()
            if index > 0 and len(tokenize) > 1:
                row = list(tokenize[1].lower())
                print(row)
                to_matrix.append(row)
    return numpy.array(to_matrix, numpy.unicode_)


start_cds_matrix = matrixer('duplexW_ZE100.aln')
stop_cds_matrix = matrixer('duplexW_EZ100.aln')
# donor_matrix = matrixer('duplexW_EI.aln')


background_state = State(DiscreteDistribution({'a': 0.25, 'c': 0.25, 'g': 0.25, 't': 0.25}), name='background')
start_codon_region_state = State(None, 'none_start_codon_state')
coding_region1 = State(DiscreteDistribution({'a': 0.25, 'c': 0.25, 'g': 0.25, 't': 0.25}), name='coding region 1')
coding_region2 = State(DiscreteDistribution({'a': 0.25, 'c': 0.25, 'g': 0.25, 't': 0.25}), name='coding region 2')
coding_region3 = State(DiscreteDistribution({'a': 0.25, 'c': 0.25, 'g': 0.25, 't': 0.25}), name='coding region 3')
stop_codon_region_state = State(None, 'none_stop_codon_state')

# donor_start = State(None, 'none_donor_start')

model = HMMWrapper()

model.add_state(background_state, 0.9)
model.add_state(start_codon_region_state)
model.add_state(coding_region1)
model.add_state(coding_region2)
model.add_state(coding_region3)
model.add_state(stop_codon_region_state)

model.add_transition(background_state, background_state, 0.9)
model.add_transition(background_state, start_codon_region_state, 0.1)
model.add_transition(coding_region1, coding_region2, 1.0)
model.add_transition(coding_region2, coding_region3, 1.0)
model.add_transition(coding_region3, coding_region1, 0.9)
model.add_transition(coding_region3, stop_codon_region_state, 0.1)

model.make_states_from_alignment(start_codon_region_state, coding_region1, start_cds_matrix, 'start codon region ')
model.make_states_from_alignment(stop_codon_region_state, background_state, stop_cds_matrix, 'stop codon region ')
# model.add_state(donor_start)
# model.add_transition(background_state, donor_start, 0.1)
# model.make_states_from_alignment(donor_start, background_state, donor_matrix, 'donor')

model.bake()

a = 'a'
c = 'c'
g = 'g'
t = 't'
seq = numpy.array([])
print(len(seq))
hmm_predictions = model.predict(seq, algorithm='viterbi')

for pre in hmm_predictions:
    print(model.states[pre].name)


def not_founder(states1, states2):
    for state1 in states1:
        found = False
        for state2 in states2:
            if state1.name == state2.name:
                found = True
        if not found:
            print(state1.name)

# print(len(model.states_before_bake))
# print(model.states)
# print('yeepi', model.model)
# not_founder(model.states_before_bake, model.states)


