import numpy
import time
from pomegranate import State
from pomegranate import DiscreteDistribution
from pomegranate import HiddenMarkovModel

example = [
    ['-', 'a', 'g', 't', 'a', 'a', 'a'],
    ['-', 'a', 'g', 't', 'a', '-', 'a'],
    ['-', 'a', 'g', 't', '-', '-', 'a'],
    ['-', 'a', 'g', 't', '-', '-', 'a'],
    ['a', 'a', 'g', 't', '-', 't', 'a'],
    ['a', '-', 'g', 't', '-', '-', 'a'],
    ['a', 'a', 'g', 't', '-', '-', 'a'],
]

data_matrix = numpy.array(example, numpy.unicode_)


class Column:
    def __init__(self, elements,):
        self.elements = elements
        self.counts = {
            'a': 0,
            'c': 0,
            'g': 0,
            't': 0,
            '-': 0,
        }
        self.delete = False
        count = (dict(zip(*numpy.unique(elements, return_counts=True))))

        for key, value in count.items():
            self.counts[key] = value

        if self.counts['-'] >= len(elements) / 2:
            self.type = 'insert'
        else:
            self.type = 'main'
            if self.counts['-']:
                self.delete = True


def column_clasify(matrix):
    classified = []
    for column in matrix.T:
        classified.append(Column(column))
    return classified


def insertion_zone(classified, index):
    zone = {
        'type': 'insert',
        'columns': []
    }
    while index < len(classified):
        if classified[index].type != 'insert':
            break
        zone['columns'].append(classified[index])
        index = index + 1
    return zone


def main_zone(data):
    group = {
        'type': 'main',
        'delete': False,
        'column': data,
    }
    if data.counts['-']:
        group['delete'] = True
    return group


def create_zones(classified):
    zones = []
    for index, data in enumerate(classified):
        if data.type == 'insert' and (not index or index and classified[index - 1].type != 'insert'):
            zones.append(insertion_zone(classified, index))
        elif data.type == 'main':
            zones.append(main_zone(data))
    return zones


def make_insert(zone, name):
    emission = {}
    total = 0
    for column in zone['columns']:
        for el in column.elements:
            if el != '-':
                if el not in emission:
                    emission[el] = 2
                    total += 2
                else:
                    emission[el] += 1
                    total += 1
    for key in emission:
        emission[key] = emission[key] / total
    # print(emission)
    return {
        'type': 'insert',
        'emission': emission,
        'zone': zone,
        'insert_state': State(DiscreteDistribution(emission), name='insert ' + name)
    }


def make_main(zone, name):
    emission = {}
    total = 0
    for el in zone['column'].elements:
        if el != '-':
            if el not in emission:
                emission[el] = 2
                total += 2
            else:
                emission[el] += 1
                total += 1

    for key in emission:
        emission[key] = emission[key] / total
    # print('main', emission)
    return {
        'type': 'main',
        'emission': emission,
        'zone': zone,
        'main_state': State(DiscreteDistribution(emission), name='main ' + name),
        'delete_state': State(None, name='none delete ' + name) if zone['delete'] else None
    }


def states(zone, name):
    if zone['type'] == 'main':
        return make_main(zone, name)
    elif zone['type'] == 'insert':
        return make_insert(zone, name)


def group_states(zones, name):
    grouped = []
    for index, subzones in enumerate(zones):
        group = states(subzones, name + str(index))
        # print(group)
        grouped.append(group)
    return grouped


def add_states(model, grouped_states):
    for group in grouped_states:
        if 'main_state' in group:
            model.add_state(group['main_state'])
        if 'delete_state' in group and group['delete_state']:
            # print('adding state', group['delete_state'].name)
            model.add_state(group['delete_state'])
        if 'insert_state' in group:
            #print(group)
            model.add_state(group['insert_state'])


def next_state(group_index, element_index, grouped_states, end_state, column_number=0):
    if group_index >= len(grouped_states):
        return end_state
    this_group = grouped_states[group_index]
    if this_group['type'] == 'main':
        element = this_group['zone']['column'].elements[element_index]
        if element != '-':
            return this_group['main_state']
        else:
            return this_group['delete_state']

    elif this_group['type'] == 'insert':
        columns = this_group['zone']['columns']
        for index, column in enumerate(columns):
            if index >= column_number:
                element = column.elements[element_index]
                if element != '-':
                    return this_group['insert_state']
        return next_state(group_index + 1, element_index, grouped_states, end_state)


def to_first(start_state, first_gstate, end_state, grouped_states):
    next_states = {
        'start': start_state,
        'total': 0,
        'states': {}
    }
    if first_gstate['type'] == 'insert':
        column = first_gstate['zone']['columns'][0]
    else:
        column = first_gstate['zone']['column']
    for index, el in enumerate(column.elements):
        next_s = next_state(0, index, grouped_states, end_state)
        if next_s.name not in next_states['states']:
            next_states['total'] += 1
            next_states['states'][next_s.name] = {
                'count': 1,
                'state': next_s,
            }
        else:
            next_states['total'] += 1
            next_states['states'][next_s.name]['count'] += 1
    return next_states


def from_main(index, end_state, grouped_states):
    init_state = grouped_states[index]
    elements = init_state['zone']['column'].elements

    main_trans = {
        'start': init_state['main_state'],
        'total': 0,
        'states': {},
    }

    delete_trans = {
        'start': init_state['delete_state'],
        'total': 0,
        'states': {},
    }

    for el_index, el in enumerate(elements):
        next_s = next_state(index + 1, el_index, grouped_states, end_state)
        if el != '-':
            main_trans['total'] += 1
            if next_s.name not in main_trans['states']:
                main_trans['total'] += 1
                main_trans['states'][next_s.name] = {
                    'count': 2,
                    'state': next_s,
                }
            else:
                main_trans['states'][next_s.name]['count'] += 1
        else:
            delete_trans['total'] += 1
            if next_s.name not in delete_trans['states']:
                delete_trans['total'] += 1
                delete_trans['states'][next_s.name] = {
                    'count': 2,
                    'state': next_s,
                }
            else:
                delete_trans['states'][next_s.name]['count'] += 1

    return main_trans, delete_trans


def from_insert(index, end_state, grouped_states):
    group = grouped_states[index]
    columns = group['zone']['columns'];
    trans = {
        'start': group['insert_state'],
        'total': 0,
        'states': {},
    }
    for col_index, column in enumerate(columns):
        for el_index, el in enumerate(column.elements):
            next_s = next_state(index, el_index, grouped_states, end_state, col_index + 1)

            if el != '-':
                trans['total'] += 1
                if next_s.name not in trans['states']:
                    trans['total'] += 1
                    trans['states'][next_s.name] = {
                        'count': 2,
                        'state': next_s,
                    }
                else:
                    trans['states'][next_s.name]['count'] += 1
    return trans


def calculate_transitions(start_state, end_state, grouped_states):
    states_from_start = to_first(start_state, grouped_states[0], end_state, grouped_states)
    all_trans = [states_from_start]
    for index, gstate in enumerate(grouped_states):
        if gstate['type'] == 'main':
            trans = from_main(index, end_state, grouped_states)
            for t in trans:
                if t['start']:
                    all_trans.append(t)
        else:
            trans = from_insert(index, end_state, grouped_states)
            all_trans.append(trans)
    return all_trans


def apply_transitions(model, transitions):
    for tran in transitions:
        for key, state_data in tran['states'].items():
            trans_prob = state_data['count'] / tran['total']
            model.add_transition(tran['start'], state_data['state'], trans_prob)


def insert_delete_main_hmm(data_matrix):
    v_columns = column_clasify(data_matrix)
    v_zones = create_zones(v_columns)
    v_grouped_states = group_states(v_zones, 'test')
    v_model = HiddenMarkovModel()
    v_first_state = State(None, name='ali_start')
    v_last_state = State(None, name='ali_end')
    v_model.add_state(v_first_state)
    v_model.add_transition(v_model.start, v_first_state, 1)
    v_model.add_state(v_last_state)
    add_states(v_model, v_grouped_states)
    v_trans = calculate_transitions(v_first_state, v_last_state, v_grouped_states)
    apply_transitions(v_model, v_trans)
    v_model.bake()
    return v_model


def insert_delete_main_hmm2(model, first_state, last_state, seq_matrix, name):
    columns = column_clasify(seq_matrix)
    zones = create_zones(columns)
    grouped_states = group_states(zones, name)
    add_states(model, grouped_states)
    trans = calculate_transitions(first_state, last_state, grouped_states)
    apply_transitions(model, trans)


class HMMWrapper:
    def __init__(self):
        self.model = HiddenMarkovModel()
        self.start = self.model.start
        self.end = self.model.end
        self.states_before_bake = []
        self.states = None

    def add_state(self, state, start_prob=0):
        self.states_before_bake.append((state, start_prob))
        self.model.add_state(state)

    def add_transition(self, start_state, end_state, prob):
        # print('adding from', start_state.name, 'to', end_state.name, prob)
        self.model.add_transition(start_state, end_state, prob)

    def bake(self):
        starter_states_no_prob = []
        free_start_prob = 1.0
        for state in self.states_before_bake:
            if 'none' not in state[0].name:
                if not state[1]:
                    starter_states_no_prob.append(state)
                else:
                    free_start_prob -= state[1]
                    print('asignado ' + str(state[1]) + ' a ' + state[0].name)
                    self.add_transition(self.start, state[0], state[1])

        len_no_prob = len(starter_states_no_prob)
        starter_prob = free_start_prob / len_no_prob
        print(len_no_prob, starter_prob)
        for state in starter_states_no_prob:
            self.add_transition(self.start, state, starter_prob)

        self.model.bake()
        self.states = self.model.states

    def make_states_from_alignment(self, first_state, last_state, seq_matrix, name):
        columns = column_clasify(seq_matrix)
        zones = create_zones(columns)
        grouped_states = group_states(zones, name)
        add_states(self, grouped_states)
        trans = calculate_transitions(first_state, last_state, grouped_states)
        apply_transitions(self, trans)

    def predict(self, *args, **kwargs):
        return self.model.predict(*args, **kwargs)
