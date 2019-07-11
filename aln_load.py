import numpy
from pomegranate import HiddenMarkovModel
from pomegranate import State
from pomegranate import DiscreteDistribution
from insertor import insert_delete_main_hmm2

filename = 'duplexW_EI.aln'
to_matrix = []
with open(filename) as file_handle:
    for index, line in enumerate(file_handle):
        tokenize = line.split()
        if index > 0 and len(tokenize) > 1:
            row = list(tokenize[1].lower())
            print(row)
            to_matrix.append(row)

data_matrix = numpy.array(to_matrix, numpy.unicode_)

background_state = State(DiscreteDistribution({'a': 0.25, 'c': 0.25, 'g': 0.25, 't': 0.25}), name='background')
donor_start = State(None, 'donor_start')

model = HiddenMarkovModel()
model.add_state(background_state)
model.add_transition(model.start, background_state, 1.0)
model.add_transition(background_state, background_state, 0.9)
model.add_transition(background_state, donor_start, 0.1)

insert_delete_main_hmm2(model, donor_start, background_state, data_matrix, 'donor')

model.bake()

print(model.states)

a = 'a'
c = 'c'
g = 'g'
t = 't'
seq = numpy.array(['a', 'g', 'g', 'g', 'c', 'c', 't', 't', 't', 'g', 'c', 't', 'g', 'g', 't', 'g', 'a', 'g', 't', 'a', 'g', 'c', 't', 'g', 't', 'c', 'c', 't', 'c', 't', 'g', 'a', 'a', 'g', 'g', 't', 'a', 'a', 'g', 'g', 'a'])
print(len(seq))
hmm_predictions = model.predict(seq, algorithm='viterbi')

for pre in hmm_predictions:
    print(model.states[pre].name)
