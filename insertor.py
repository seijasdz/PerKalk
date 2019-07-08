import numpy
from pomegranate import State
from pomegranate import DiscreteDistribution
from pomegranate import HiddenMarkovModel

example = [
    ['-', 'a', 'a', '-', 'a'],
    ['-', 'a', 'a', 'a', 'a'],
    ['-', '-', 'a', '-', '-'],
    ['-', 'a', 'a', 'a', '-'],
    ['a', '-', 'a', '-', '-'],
    ['-', 'a', 'a', '-', '-'],
    ['a', 'a', 'a', '-', '-'],
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
    emission = {
        'a': 1,
        'c': 1,
        'g': 1,
        't': 1,
    }
    total = 4
    for column in zone['columns']:
        for el in column.elements:
            if el != '-':
                emission[el] = emission[el] + 1
                total = total + 1
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
    emission = {
        'a': 1,
        'c': 1,
        'g': 1,
        't': 1,
    }
    total = 4
    for el in zone['column'].elements:
        if el != '-':
            emission[el] = emission[el] + 1
            total = total + 1

    for key in emission:
        emission[key] = emission[key] / total
    # print('main', emission)
    return {
        'type': 'main',
        'emission': emission,
        'zone': zone,
        'main_state': State(DiscreteDistribution(emission), name='main ' + name),
        'delete_state': State(None, name='delete ' + name) if zone['delete'] else None
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
            model.add_state(group['delete_state'])
        if 'insert_state' in group:
            model.add_state(group['insert_state'])


def next_state(group_index, element_index, grouped_states, end_state):
    this_group = grouped_states[group_index]
    if this_group['type'] == 'main':
        element = this_group['zone']['column'].elements[element_index]
        if element != '-':
            return this_group['main_state']
        else:
            return this_group['delete_state']

    elif this_group['type'] == 'insert':
        columns = this_group['zone']['columns']
        for column in columns:
            element = column.elements[element_index]
            if element != '-':
                return this_group['insert_state']
        if group_index + 1 < len(grouped_states):
            return next_state(group_index + 1, element_index, grouped_states, end_state)
        return end_state


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
                # 'state': next_s,
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

    for el_index, el in enumerate(elements):
        next_s = next_state(index + 1, el_index, grouped_states, end_state)
        if el != '-':
            main_trans['total'] += 1
            if next_s.name not in main_trans['states']:
                main_trans['total'] += 1
                main_trans['states'][next_s.name] = {
                    'count': 2,
                    # 'state:' next_s,
                }
            else:
                main_trans['states'][next_s.name]['count'] += 1
        else:
            print('coming soon')

    print(main_trans)





def transitions(model, start_state, end_state, grouped_states):
    states_from_start = to_first(start_state, grouped_states[0], end_state, grouped_states)
    for index, gstate in enumerate(grouped_states):
        if gstate['type'] == 'main':
            transitions = from_main(index, end_state, grouped_states)
        else:
            print('i')

v_columns = column_clasify(data_matrix)

v_zones = create_zones(v_columns)

v_grouped_states = group_states(v_zones, 'test')

v_model = HiddenMarkovModel()
v_first_state = State(None, name='ali_start')
v_last_state = State(None, name='ali_end')

v_model.add_state(v_first_state)
v_model.add_state(v_last_state)

add_states(v_model, v_grouped_states)
transitions(v_model,v_first_state, v_last_state, v_grouped_states)


