from pomegranate import DiscreteDistribution
from pomegranate import State
from converter_to import converter_to

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


def sequence_state_factory(states_data, name):
    states = []
    for index, data in enumerate(states_data):
        state = State(DiscreteDistribution(data.states_distribution), name=name + str(index))
        states.append(state)
    return states


def spacer_states_maker(quantity, distribution, name):
    states = []
    for i in range(0, quantity):
        state = State(DiscreteDistribution(distribution), name=name + str(i))
        states.append(state)
    return states


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


def add_sequence(model, states):
    for state in states:
        model.add_state(state)

    for index, state in enumerate(states):
        if index < len(states) - 1:
            model.add_transition(state, states[index + 1], 1.0)


def add_variable_length_sequence(model, states, end_state):
    for state in states:
        model.add_state(state)

    for index, state in enumerate(states):
        if index < len(states) - 1:
            model.add_transition(state, states[index + 1], 0.8)
            model.add_transition(state, end_state, 0.2)
    model.add_transition(states[-1], end_state, 1.0)


def load_long_training_examples(filename, n):
    with open(filename) as file:
        first = file.read().replace('\n', '').replace(' ', '').split('>')[1:n]
        second = [converter_to(f.split('.')[1].lower(), 2) for f in first]
    return second


class StateNotFoundException(Exception):
    pass


def get_state(model, name):
    for state in model.states:
        if state.name == name:
            return state
    raise StateNotFoundException('State not found ' + name)


if __name__ == '__main__':
    load_long_training_examples('all_promoters_-499_5.fa', 4002)
