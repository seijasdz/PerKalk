import numpy


example = [
    ['-', 'a', 'a', '-', 'a'],
    ['-', 'a', 'a', 'a', '-'],
    ['-', '-', 'a', '-', '-'],
    ['-', 'a', 'a', 'a', '-'],
    ['a', '-', 'a', '-', '-'],
    ['-', 'a', 'a', '-', '-'],
    ['a', 'a', 'a', '-', 'a'],
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
    }
    m_zone = {
        'type': 'main',
        'columns': [data]
    }
    group['main'] = m_zone
    if data.counts['-']:
        e_zone = {
            'type': 'elimination',
            'columns': [data]
        }
        group['elimination'] = e_zone
    return group


def create_zones(classified):
    zones = []
    for index, data in enumerate(classified):
        if data.type == 'insert' and (not index or index and classified[index - 1].type != 'insert'):
            zones.append(insertion_zone(classified, index))
        elif data.type == 'main':
            zones.append(main_zone(data))
    return zones


def make_insert(zone):
    print('haciendo insert')
    print(zone)
    return {}


def make_main(zone):
    print('haciendo main')
    return {}


def states(zone):
    if zone['type'] == 'main':
        return make_main(zone)
    elif zone['type'] == 'insert':
        return make_insert(zone)


def group_states(zones):
    grouped = []
    for index, subzones in enumerate(zones):
        group = states(subzones)
        grouped.append(group)
    return grouped


v_columns = column_clasify(data_matrix)
v_zones = create_zones(v_columns)
v_grouped_states = group_states(v_zones)
