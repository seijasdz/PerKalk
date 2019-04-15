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

    def add_counts(self, nucleotides):
        for nuc in nucleotides:
            self.add_count(nuc)

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

        def cond_add_state(id, dist, name):
            if self.model:
                self.states[id] = self.model.add_state(dist, name)

        self.emission_probabilities = {
            'A': self.counts['A'] / self.nucleotides_total if self.nucleotides_total else 0,
            'C': self.counts['C'] / self.nucleotides_total if self.nucleotides_total else 0,
            'G': self.counts['G'] / self.nucleotides_total if self.nucleotides_total else 0,
            'T': self.counts['T'] / self.nucleotides_total if self.nucleotides_total else 0,
        }
        print(self.emission_probabilities)

        if self.type == StateFactory.INSERT:

            cond_add_state('insert', DiscreteDistribution(self.emission_probabilities),
                                                         'i_st_' + str(self.idx))
            for tr in self.trans:
                if tr[4] == StateFactory.INSERT:
                    insert_to_self = insert_to_self + 1
                elif tr[4] == StateFactory.MAIN:
                    insert_to_main = insert_to_main + 1
                elif tr[4] == StateFactory.DELETE:
                    insert_to_delete = insert_to_delete + 1

        elif self.type == StateFactory.MAIN:
            cond_add_state('main', DiscreteDistribution(self.emission_probabilities),
                                                         'm_st_' + str(self.idx))
            if self.counts['-']:
                cond_add_state('delete', None, 'd_st_' + str(self.idx))
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
            'main_to_insert': main_to_insert / total_main_trans if self.type == StateFactory.MAIN and total_main_trans else 0,
            'main_to_main': main_to_main / total_main_trans if self.type == StateFactory.MAIN and total_main_trans else 0,
            'main_to_del': main_to_del / total_main_trans if self.type == StateFactory.MAIN and total_main_trans else 0,
            'del_to_insert': del_to_insert / total_delete_trans if total_delete_trans else 0,
            'del_to_main': del_to_main / total_delete_trans if total_delete_trans else 0,
            'del_to_del': del_to_del / total_delete_trans if total_delete_trans else 0,
        }

        print(self.transition_probabilities)


class StateFactory:

    INSERT = 0
    MAIN = 1
    DELETE = 2

    def __init__(self, model, previous_states, next_states):
        self.columns = []
        self.zones = []
        self.model = model
        self.previous_states = previous_states
        self.next_states = next_states

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

        first_zone = Zone(StateFactory.MAIN, None, None)
        first_zone.add_counts(['A', 'A', 'A'])

        for row, element in enumerate(self.columns[0]['info']['content']):
            trans = check_next(row, 0)
            first_zone.add_transition((row, trans[2], element, trans[1], trans[0]))

        self.zones.append(first_zone)


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

    def make_weave_states(self):
        for idx, zone in enumerate(self.zones):
                print('in zone', idx)
                for key, state in zone.states.items():
                    print(key, state)




lines = []
file_name = 'test'
aln_file_name = file_name + '.aln'

with open(aln_file_name) as aln_file:
    for line in aln_file:
        if 'CLUSTAL' not in line and len(line) > 1 and line[0] != ' ':
            thisline = list(line[16:-1])
            lines.append(thisline)
            print(thisline)

array = numpy.array(lines, numpy.unicode_)
ref = array[:, 0]

rows = array.shape[0]
print('rows', rows)

def same_dist():
    return DiscreteDistribution({'a': 0.25, 'c': 0.25, 'g': 0.25, 't': 0.25})


counts = []
model = ModelWrapper()
background_state = model.add_state(same_dist(), 'bk')

model.model.add_transition(model.model.start, background_state, 0.5)
model.model.add_transition(background_state, background_state, 0.9)

st_fac = StateFactory(model, [(background_state, 0.1), (model.model.start, 0.5)], [(background_state, 1)])

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
st_fac.make_weave_states()
