import numpy
from pomegranate import State
from pomegranate import DiscreteDistribution
from insertor import HMMWrapper
from gene_sample_extractor import seqs_from

def matrixer(filename):
    to_matrix = []
    with open(filename) as file_handle:
        print(filename)
        for index, line in enumerate(file_handle):
            tokenize = line.split()
            if index > 0 and len(tokenize) > 1:
                row = list(tokenize[1].lower())
                print(row)
                to_matrix.append(row)
    return numpy.array(to_matrix, numpy.unicode_)


start_cds_matrix = matrixer('duplexW_ZE100.aln')
stop_cds_matrix = matrixer('duplexW_EZ100.aln')

donor0_matrix = matrixer('donor0.aln')
donor1_matrix = matrixer('donor1.aln')
donor2_matrix = matrixer('donor2.aln')

acceptor0_matrix = matrixer('acceptor0.aln')
acceptor1_matrix = matrixer('acceptor1.aln')
acceptor2_matrix = matrixer('acceptor2.aln')


# donor_matrix = matrixer('duplexW_EI.aln')


background_state = State(DiscreteDistribution({'a': 0.25, 'c': 0.25, 'g': 0.25, 't': 0.25}), name='background')
start_codon_region_state = State(None, name='none_start_codon_state')
coding_region1 = State(DiscreteDistribution({'a': 0.25, 'c': 0.25, 'g': 0.25, 't': 0.25}), name='coding region 1')
coding_region2 = State(DiscreteDistribution({'a': 0.25, 'c': 0.25, 'g': 0.25, 't': 0.25}), name='coding region 2')
coding_region3 = State(DiscreteDistribution({'a': 0.25, 'c': 0.25, 'g': 0.25, 't': 0.25}), name='coding region 3')
stop_codon_region_state = State(None, name='none_stop_codon_state')
donor_start = State(None, name='none donor region')
donor_type0_start = State(None, name='none donor0 region start')
donor_type1_start = State(None, name='none donor1 region start')
donor_type2_start = State(None, name='none donor2 region start')
interior_intron0 = State(DiscreteDistribution({'a': 0.25, 'c': 0.25, 'g': 0.25, 't': 0.25}), name='interior_intron0')
interior_intron1 = State(DiscreteDistribution({'a': 0.25, 'c': 0.25, 'g': 0.25, 't': 0.25}), name='interior_intron1')
interior_intron2 = State(DiscreteDistribution({'a': 0.25, 'c': 0.25, 'g': 0.25, 't': 0.25}), name='interior_intron2')
acceptor0_start = State(None, name='none acceptor0 region start')
acceptor1_start = State(None, name='none acceptor1 region start')
acceptor2_start = State(None, name='none acceptor2 region start')


model = HMMWrapper()

model.add_state(background_state, 0.9)
model.add_state(start_codon_region_state)
model.add_state(coding_region1)
model.add_state(coding_region2)
model.add_state(coding_region3)
model.add_state(stop_codon_region_state)
model.add_state(donor_start)
model.add_state(donor_type0_start)
model.add_state(donor_type1_start)
model.add_state(donor_type2_start)
model.add_state(interior_intron0)
model.add_state(interior_intron1)
model.add_state(interior_intron2)
model.add_state(acceptor0_start)
model.add_state(acceptor1_start)
model.add_state(acceptor2_start)


model.add_transition(background_state, background_state, 0.99)
model.add_transition(background_state, start_codon_region_state, 0.01)
model.add_transition(coding_region1, coding_region2, 1.0)
model.add_transition(coding_region2, coding_region3, 1.0)
model.add_transition(coding_region3, coding_region1, 0.8)
model.add_transition(coding_region3, stop_codon_region_state, 0.01)
model.add_transition(coding_region3, donor_start, 0.19)
model.add_transition(donor_start, donor_type0_start, 0.34)
model.add_transition(donor_start, donor_type1_start, 0.33)
model.add_transition(donor_start, donor_type2_start, 0.33)
model.add_transition(interior_intron0, interior_intron0, 0.95)
model.add_transition(interior_intron0, acceptor0_start, 0.05)
model.add_transition(interior_intron1, interior_intron1, 0.95)
model.add_transition(interior_intron1, acceptor1_start, 0.05)
model.add_transition(interior_intron2, interior_intron2, 0.95)
model.add_transition(interior_intron2, acceptor2_start, 0.05)


model.make_states_from_alignment(start_codon_region_state, coding_region1, start_cds_matrix, 'start codon region ')
model.make_states_from_alignment(stop_codon_region_state, background_state, stop_cds_matrix, 'stop codon region ')

model.make_states_from_alignment(donor_type0_start, interior_intron0, donor0_matrix, 'donor type0 region')
model.make_states_from_alignment(acceptor0_start, coding_region1, acceptor0_matrix, 'acceptor type0 region')

model.make_states_from_alignment(donor_type1_start, interior_intron1, donor1_matrix, 'donor type1 region')
model.make_states_from_alignment(acceptor1_start, coding_region1, acceptor1_matrix, 'acceptor type1 region')

model.make_states_from_alignment(donor_type2_start, interior_intron2, donor2_matrix, 'donor type2 region')
model.make_states_from_alignment(acceptor2_start, coding_region1, acceptor2_matrix, 'acceptor type2 region')

model.bake()

a = 'a'
c = 'c'
g = 'g'
t = 't'
s = seqs_from('example_seq.ebi')[1]

seqx = list(s[36000: 46000])
print(seqx)
print(len(s))
seq = numpy.array(seqx)
print(len(seq))
hmm_predictions = model.predict(seq, algorithm='viterbi')

count = 0
for pre in hmm_predictions:
    if 'none' not in model.states[pre].name and 'None' not in model.states[pre].name and 'back' not in model.states[pre].name:
        print(model.states[pre].name)
        print(count, seq[count])
        count += 1


def not_founder(states1, states2):
    for state1 in states1:
        found = False
        for state2 in states2:
            if state1.name == state2.name:
                found = True
        if not found:
            print(state1.name)


with open('out.txt', 'w') as out_file:
    out_file.write(model.model.to_json())

# print(len(model.states_before_bake))
# print(model.states)
# print('yeepi', model.model)
# not_founder(model.states_before_bake, model.states)


