import numpy


example = [
    ['-', 'a', 'a', 'a', 'a'],
    ['-', 'a', 'a', 'a', 'a'],
    ['-', '-', 'a', 'a', 'a'],
    ['-', 'a', 'a', 'a', 'a'],
    ['a', '-', 'a', 'a', 'a'],
    ['-', 'a', 'a', 'a', 'a'],
    ['a', 'a', 'a', 'a', 'a'],
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


def create_zones(classified):
    zones = []
    for index, column in enumerate(classified):




classified_columns = column_clasify(data_matrix)
create_zones(classified_columns)
