from matrix_from_aln import create_matrix
from matrix_from_aln import load_to_fit
from insertor import HMMWrapper
from pomegranate import State
from pomegranate import DiscreteDistribution
import numpy

tss_matrix = create_matrix('TSS.aln')

background_state = State(DiscreteDistribution({'a': 0.25, 'c': 0.25, 'g': 0.25, 't': 0.25}), name='background')
start_tss = State(None, name='none start tss')

model = HMMWrapper()

model.add_state(background_state, 0.9)
model.add_state(start_tss, 0.1)

model.add_transition(model.start, background_state, 1.0)
model.add_transition(background_state, background_state, 0.9)
model.add_transition(background_state, start_tss, 0.1)

model.make_states_from_alignment(start_tss, background_state, tss_matrix, 'tss de prueba')

model.bake()

a = 'a'
c = 'c'
g = 'g'
t = 't'
seq = numpy.array(['a', 'g', 't', 'c', 't', 'c', 'a', 'c', 't', 'c', 'c', 'g', 't', 'c', 't', 'g', 'a', 'a', 't', 't', 'c', 'c', 't', 'c', 't', 'c', 'c', 'g', 't', 'c', 't', 'c', 'c', 'c', 't', 'c', 'c', 'c', 'a', 'c', 'c', 'c', 'c', 'g', 'g', 'c', 'c', 'g', 't', 'c', 't', 'a', 't', 'g', 'c', 't', 'c', 'c', 'a', 'g', 'g', 'c', 'c', 'c', 't', 'c', 't', 'c', 'c', 't', 'c', 'g', 'c', 'g', 'g'])
print(len(seq))
hmm_predictions = model.predict(seq, algorithm='viterbi')

for pre in hmm_predictions:
    print(model.states[pre].name)

model.model.fit(load_to_fit('TSS.f'))


print(len(seq))
hmm_predictions = model.predict(seq, algorithm='viterbi')

for pre in hmm_predictions:
    print(model.states[pre].name)

