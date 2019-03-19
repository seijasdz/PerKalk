from pomegranate import State
from pomegranate import DiscreteDistribution
from pomegranate import HiddenMarkovModel


class ModelWrapper:
    def __init__(self):
        self.model = HiddenMarkovModel()

    def add_state(self, distribution, name):
        state = State(distribution, name=name)
        self.model.add_state(state)
        return state

    def bake(self):
        self.model.bake()

    def viterbi(self, seq):
        return self.model.viterbi(seq)

    def add_transition(self, states, next_state_data):
        for state in states:
            for next_data in next_state_data:
                self.model.add_transition(state, next_data[0], next_data[1])


def with_variations(dist, name):
    st = State(dist, name=name)
    sti = State(DiscreteDistribution({'a': 0.25, 'c': 0.25, 'g': 0.25, 't': 0.25}), name='i_' + name)
    std = State(None, name='d_' + name)
    return st, sti, std


def transition_maker(hmmodel, states, next_state_data):
    for state in states:
        for next_data in next_state_data:
            hmmodel.add_transition(state, next_data[0], next_data[1])
