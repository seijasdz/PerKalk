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
    return [zone]

def main_zone(data):
    zones = []
    m_zone = {
        'type': 'main',
        'columns': [data]
    }
    zones.append(m_zone)
    if data.counts['-']:
        print('eli')
        e_zone = {
            'type': 'elimination',
            'columns': [data]
        }
        zones.append(e_zone)
    return zones


def create_zones(classified):
    zones = []
    for index, data in enumerate(classified):
        print(data.type)
        if data.type == 'insert' and (not index or index and classified[index - 1].type != 'insert'):
            zones.append(insertion_zone(classified, index))
        elif data.type == 'main':
            zones.append(main_zone(data))
    print(zones)


classified_columns = column_clasify(data_matrix)
create_zones(classified_columns)
