import numpy
from pomegranate import State
from pomegranate import DiscreteDistribution
from pomegranate import HiddenMarkovModel
from converter_to_two import converter_to
from matrix_from_aln import matrix_from_fasta
from gene_sample_extractor import seqs_from
import itertools
import exon_calculator


example = [
    ['a', 'c', 't', 'a', 't', 'g', 'c', 'a', 'a'],
    ['a', 'c', 't', 'a', 't', 'g', 'c', 'a', 'a'],
    ['a', 'c', 't', 'a', 't', 'g', 'c', 'a', 'a'],
    ['a', 'c', 't', 'a', 't', 'g', 'c', 'a', 'a'],
    ['a', 'c', 't', 'a', 't', 'g', 'c', 'a', 'a'],
    ['a', 'c', 't', 'a', 't', 'g', 'c', 'a', 'a'],
    ['a', 'c', 't', 'a', 't', 'g', 'c', 'a', 'a'],
    ['a', 'c', 't', 'a', 't', 'g', 'c', 'a', 'a'],
    ['t', 'c', 'g', 'a', 't', 'g', 'c', 'a', 'a'],
    ['t', 'c', 'g', 'a', 't', 'g', 'c', 'a', 'a'],
    ['t', 'c', 'g', 'a', 't', 'g', 'c', 'a', 'a'],
    ['t', 'c', 'g', 'a', 't', 'g', 'c', 'a', 'a'],
    ['t', 'c', 'g', 'g', 't', 'g', 'c', 'a', 'a'],
    ['t', 'c', 'g', 'a', 't', 'g', 'c', 'a', 'a'],
    ['t', 'c', 'g', 'a', 't', 'g', 'c', 'a', 'a'],

]


class HighOrderState:
    def __init__(self, columns):
        self.states_distribution = self.calculate_states(columns)

    @staticmethod
    def calculate_states(columns):
        states_distribution = {}
        for index, base1 in enumerate(columns[0]):
            name = base1 + '|'
            for i, col in enumerate(columns):
                if i:
                    name += col[index]
            if name not in states_distribution:
                states_distribution[name] = 1
            else:
                states_distribution[name] += 1

        for key, value in states_distribution.items():
            states_distribution[key] /= len(columns[0])
        return states_distribution


def classify(matrix, order=1):
    state_data = []
    for index, column in enumerate(matrix.T):
        columns = [column]
        for x in range(order, 0, -1):
            col = matrix.T[index - x]
            columns.append(col)
        if index > order - 1:
            data = HighOrderState(columns)
            state_data.append(data)
    return state_data


def sequence_state_factory(states_data, name):
    states = []
    for index, data in enumerate(states_data):
        state = State(DiscreteDistribution(data.states_distribution), name=name + str(index))
        states.append(state)
    return states


def add_sequence(model, states):
    for state in states:
        model.add_state(state)

    for index, state in enumerate(states):
        if index < len(states) - 1:
            model.add_transition(state, states[index + 1], 1.0)


def foo(l):
    yield from itertools.product(*([l] * 3))


def temporal_back(order):
    states = {}
    back_prob = 1.0 / pow(4, order + 1)
    nucleotides = ['a', 'c', 'g', 't']
    t = foo(nucleotides)
    for s in t:
        state = ''.join(s)
        states[state[:1] + '|' + state[1:]] = back_prob
    return states


v_lines = exon_calculator.get_corrected_lines('exonesW.txt')
c0, c1, c2 = exon_calculator.calculate_proba(v_lines)




data_matrix = numpy.array(matrix_from_fasta('duplexW_ZE100.f'))
# data_matrix = numpy.array(example, numpy.unicode)

v_states_data = classify(data_matrix, 2)
v_states = sequence_state_factory(v_states_data, 'start zone')

model = HiddenMarkovModel()

back = State(DiscreteDistribution(temporal_back(2)), name='back')
coding_state0 = State(DiscreteDistribution(c0.p), 'coding state 0')
coding_state1 = State(DiscreteDistribution(c1.p), 'coding state 1')
coding_state2 = State(DiscreteDistribution(c2.p), 'coding state 2')

model.add_state(back)
model.add_state(coding_state0)
model.add_state(coding_state1)
model.add_state(coding_state2)


model.add_transition(model.start, back, 1.0)
model.add_transition(back, back, 0.999)

model.add_transition(coding_state0, coding_state1, 1.0)
model.add_transition(coding_state1, coding_state2, 1.0)
print(c2.trans_probs['default'], c2.trans_probs['end'])
model.add_transition(coding_state2, coding_state0, c2.trans_probs['default'])
model.add_transition(coding_state2, back, c2.trans_probs['end'])


add_sequence(model, v_states)

model.add_transition(back, v_states[0], 0.001)
model.add_transition(v_states[-1], coding_state0, 1.0)

model.bake()


string = seqs_from('sequence_body.ebi')[0][0:100000]
print(string)
test_seq = list(string)
two_seq = converter_to(test_seq, 2)
print(two_seq)
seq = numpy.array(two_seq, numpy.unicode_)

hmm_predictions = model.predict(seq, algorithm='viterbi')


count = 0
for i, pre in enumerate(hmm_predictions):
    if 'none' not in model.states[pre].name and 'None' not in model.states[pre].name and 'back' not in model.states[pre].name:
        print(model.states[pre].name, seq[i - 1])
        #print(count, seq[count])
        count += 1

