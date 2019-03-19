import numpy
from helpers import ModelWrapper
from pomegranate import DiscreteDistribution

chars = ['-', 'A', 'C', 'G', 'T']


class Zone:
    def __init__(self, type, model, idx):
        self.type = type
        self.trans = set()
        self.counts = {
            'A': 0,
            'C': 0,
            'G': 0,
            'T': 0,
            '-': 0,
        }
        self.idx = idx
        self.model = model
        self.states = {}
        self.transition_probabilities = {}
        self.emission_probabilities = {}
        self.nucleotides_total = 0
        self.states = {}

    def add_count(self, nucleotide):
        if nucleotide != '-':
            self.nucleotides_total = self.nucleotides_total + 1
        self.counts[nucleotide] = self.counts[nucleotide] + 1

    def add_transition(self, trans):
        self.trans.add(trans)

    def calc_trans_probs(self):
        insert_to_self = 0
        insert_to_main = 0
        insert_to_delete = 0
        main_to_insert = 0
        main_to_main = 0
        main_to_del = 0
        del_to_insert = 0
        del_to_main = 0
        del_to_del = 0
        total_trans = len(self.trans)
        total_main_trans = 0
        total_delete_trans = 0

        self.emission_probabilities = {
            'A': self.counts['A'] / self.nucleotides_total,
            'C': self.counts['C'] / self.nucleotides_total,
            'G': self.counts['G'] / self.nucleotides_total,
            'T': self.counts['T'] / self.nucleotides_total,
        }
        print(self.emission_probabilities)

        if self.type == StateFactory.INSERT:
            self.states['insert'] = self.model.add_state(DiscreteDistribution(self.emission_probabilities),
                                                         'i_st_' + str(self.idx))
            for tr in self.trans:
                if tr[4] == StateFactory.INSERT:
                    insert_to_self = insert_to_self + 1
                elif tr[4] == StateFactory.MAIN:
                    insert_to_main = insert_to_main + 1
                elif tr[4] == StateFactory.DELETE:
                    insert_to_delete = insert_to_delete + 1

        elif self.type == StateFactory.MAIN:
            self.states['main'] = self.model.add_state(DiscreteDistribution(self.emission_probabilities),
                                                         'm_st_' + str(self.idx))
            if self.counts['-']:
                print('agregado delete')
                self.states['main'] = self.model.add_state(None, 'd_st_' + str(self.idx))
            for tr in self.trans:
                if tr[4] == StateFactory.INSERT:
                    if tr[2] == '-':
                        del_to_insert = del_to_insert + 1
                    else:
                        main_to_insert = main_to_insert + 1
                if tr[4] == StateFactory.MAIN:
                    if tr[2] == '-':
                        del_to_main = del_to_main + 1
                    else:
                        main_to_main = main_to_main + 1
                if tr[4] == StateFactory.DELETE:
                    if tr[2] == '-':
                        del_to_del = del_to_del + 1
                    else:
                        main_to_del = main_to_del + 1

            for tran in self.trans:
                if tran[2] == '-':
                    total_delete_trans = total_delete_trans + 1
                else:
                    total_main_trans = total_main_trans + 1

        self.transition_probabilities = {
            'insert_to_self': insert_to_self / total_trans if self.type == StateFactory.INSERT else 0,
            'insert_to_main': insert_to_main / total_trans if self.type == StateFactory.INSERT else 0,
            'insert_to_delete': insert_to_delete / total_trans if self.type == StateFactory.INSERT else 0,
            'main_to_insert': main_to_insert / total_main_trans if self.type == StateFactory.MAIN else 0,
            'main_to_main': main_to_main / total_main_trans if self.type == StateFactory.MAIN else 0,
            'main_to_del': main_to_del / total_main_trans if self.type == StateFactory.MAIN else 0,
            'del_to_insert': del_to_insert / total_delete_trans if total_delete_trans else 0,
            'del_to_main': del_to_main / total_delete_trans if total_delete_trans else 0,
            'del_to_del': del_to_del / total_delete_trans if total_delete_trans else 0,
        }

        print(self.transition_probabilities)


class StateFactory:

    INSERT = 0
    MAIN = 1
    DELETE = 2

    def __init__(self, model, prev, next):
        self.columns = []
        self.zones = []
        self.model = model

    def add_info(self, info):
        if info['count']['-'] >= len(info['content']) / 2:
            column = {
                'type': StateFactory.INSERT,
                'info': info
            }
            self.columns.append(column)
        else:
            column = {
                'type': StateFactory.MAIN,
                'info': info
            }
            self.columns.append(column)

    def make_zones(self):
        def check_next(row, col):
            if col >= len(self.columns):
                return StateFactory.MAIN, 'X', col
            if self.columns[col]['type'] == StateFactory.MAIN:
                if self.columns[col]['info']['content'][row] != '-':
                    return self.columns[col]['type'], self.columns[col]['info']['content'][row], col
                else:
                    return StateFactory.DELETE, '-', col
            elif self.columns[col]['type'] == StateFactory.INSERT:
                if self.columns[col]['info']['content'][row] == '-':
                    return check_next(row, col + 1)
                else:
                    return self.columns[col]['type'], self.columns[col]['info']['content'][row], col

        for idx, col in enumerate(self.columns):
            if col['type'] == StateFactory.INSERT:
                if not idx or idx and self.columns[idx - 1]['type'] != StateFactory.INSERT:
                    self.zones.append(Zone(StateFactory.INSERT, self.model, idx))

                zone = self.zones[-1]
                for row, element in enumerate(col['info']['content']):
                    zone.add_count(element)
                    trans = check_next(row, idx + 1)
                    if element != '-':
                        zone.add_transition((row, trans[2], element, trans[1], trans[0]))

                print(zone.trans, zone.counts)
            elif col['type'] == StateFactory.MAIN:
                self.zones.append(Zone(StateFactory.MAIN, self.model, idx))
                zone = self.zones[-1]
                for row_idx, element in enumerate(col['info']['content']):
                    zone.add_count(element)
                    trans = check_next(row_idx, idx + 1)
                    zone.add_transition((row_idx, trans[2], element, trans[1], trans[0]))
                print(zone.trans, zone.counts)

        for zone in self.zones:
            zone.calc_trans_probs()




lines = []
file_name = 'test'
aln_file_name = file_name + '.aln'

with open(aln_file_name) as aln_file:
    for line in aln_file:
        if 'CLUSTAL' not in line and len(line) > 1 and line[0] != ' ':
            lines.append(list(line[16:-1]))
            print(list(line[16:-1]))

array = numpy.array(lines, numpy.unicode_)
ref = array[:, 0]

rows = array.shape[0]
print('rows', rows)

counts = []
model = ModelWrapper()
st_fac = StateFactory(model, None, None)

for col_content in array.T:
    count = (dict(zip(*numpy.unique(col_content, return_counts=True))))
    for char in chars:
        if char not in count:
            count[char] = 0
    info = {
        'content': col_content,
        'count': count
    }
    st_fac.add_info(info)

st_fac.make_zones()

