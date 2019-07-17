import numpy
from pomegranate import State
from pomegranate import DiscreteDistribution
from pomegranate import HiddenMarkovModel
from converter_to_two import converter_to_two


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

data_matrix = numpy.array(example, numpy.unicode_)


class FirstOrderState:
    def __init__(self, column1, column0):
        self.states_distribution = self.calculate_states(column1, column0)

    @staticmethod
    def calculate_states(column1, column0):
        states_distribution = {}
        for index, base1 in enumerate(column1):
            base0 = column0[index]
            name = base1 + str('|') + base0
            if name not in states_distribution:
                states_distribution[name] = 1
            else:
                states_distribution[name] += 1

        for key, value in states_distribution.items():
            states_distribution[key] /= len(column1)
        return states_distribution


def classify(matrix):
    state_data = []
    for index, column in enumerate(matrix.T):
        if index > 0:
            data = FirstOrderState(column, matrix.T[index - 1])
            state_data.append(data)
    return state_data


def sequence_state_factory(states_data, name):
    states = []
    for index, data in enumerate(states_data):
        print(data.states_distribution)
        state = State(DiscreteDistribution(data.states_distribution), name=name + str(index))
        states.append(state)
    return states


def add_sequence(model, states):
    for state in states:
        model.add_state(state)

    for index, state in enumerate(states):
        if index < len(states) - 1:
            print('trans from', state.name, 'to', states[index + 1].name)
            model.add_transition(state, states[index + 1], 1.0)


v_states_data = classify(data_matrix)
v_states = sequence_state_factory(v_states_data, 'test')


model = HiddenMarkovModel()

back_prob = 1 / 16
back_states = {
    'a|a': back_prob,
    'a|c': back_prob,
    'a|g': back_prob,
    'a|t': back_prob,
    'c|a': back_prob,
    'c|c': back_prob,
    'c|g': back_prob,
    'c|t': back_prob,
    'g|a': back_prob,
    'g|c': back_prob,
    'g|g': back_prob,
    'g|t': back_prob,
    't|a': back_prob,
    't|c': back_prob,
    't|g': back_prob,
    't|t': back_prob,
}

back = State(DiscreteDistribution(back_states), name='back')
model.add_state(back)
model.add_transition(model.start, back, 1.0)
model.add_transition(back, back, 0.1)

add_sequence(model, v_states)

model.add_transition(back, v_states[0], 0.9)
model.add_transition(v_states[-1], back, 1.0)

model.bake()

print(model.states)

test_seq = list('acgtacgtacgtacgtacgtacgtacgtactatgcaa')
two_seq = converter_to_two(test_seq)
print(two_seq)
seq = numpy.array(two_seq, numpy.unicode_)
print(seq)
hmm_predictions = model.predict(seq, algorithm='viterbi')

#print(model.states)

count = 0
for pre in hmm_predictions:
    #if 'none' not in model.states[pre].name and 'None' not in model.states[pre].name and 'back' not in model.states[pre].name:
        print(model.states[pre].name)
        #print(count, seq[count])
        #count += 1