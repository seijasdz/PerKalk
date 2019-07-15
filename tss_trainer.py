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

model.add_state(background_state)
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
seq = numpy.array(['g', 'a', 't', 'c', 'a', 't', 'g', 'c', 'c', 'c', 'a', 't', 'c', 't', 'c', 'a', 'c', 'g', 'c', 'c', 'a', 'g', 'g', 'g', 'g', 'g', 'c', 'g', 'g', 'a'])
print(len(seq))
hmm_predictions = model.predict(seq, algorithm='viterbi')

for pre in hmm_predictions:
    print(model.states[pre].name)

model.model.fit(load_to_fit('TSS.f'))


print(len(seq))
hmm_predictions = model.predict(seq, algorithm='viterbi')

logp, path = model.model.viterbi(seq)

print(logp, path)

for pre in hmm_predictions:
    print(model.states[pre].name)

