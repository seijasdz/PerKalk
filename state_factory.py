import numpy


chars = ['-', 'A', 'C', 'G', 'T']

class StateFactory:

    INSERT_STATE = 0
    MAIN_STATE = 1
    DELETE_STATE = 2

    def __init__(self, prev, next):
        self.columns = []
        self.zones = []

    def add_info(self, info):
        if info['count']['-'] >= len(info['content']) / 2:
            column = {
                'type': StateFactory.INSERT_STATE,
                'info': info
            }
            self.columns.append(column)
        else:
            column = {
                'type': StateFactory.MAIN_STATE,
                'info': info
            }
            self.columns.append(column)

    def make_zones(self):
        def check_next(row, col):
            if col >= len(self.columns):
                return 'X'
            if self.columns[col]['type'] == StateFactory.MAIN_STATE:
                if self.columns[col]['info']['content'][row] != '-':
                    return self.columns[col]['type'], self.columns[col]['info']['content'][row]
                else:
                    return StateFactory.DELETE_STATE, '-'
            elif self.columns[col]['type'] == StateFactory.INSERT_STATE:
                if self.columns[col]['info']['content'][row] == '-':
                    return check_next(row, col + 1)
                else:
                    return self.columns[col]['type'], self.columns[col]['info']['content'][row]

        for idx, col in enumerate(self.columns):
            if col['type'] == StateFactory.INSERT_STATE and\
                    (not idx or idx and self.columns[idx - 1]['type'] != StateFactory.INSERT_STATE):
                for row, elem in enumerate(col['info']['content']):
                    print(col, check_next(row, idx + 1))




    def makeold(self):
        for index, info in enumerate(self.infos):
            if info['type'] == StateFactory.INSERT_STATE and\
                    (not index or index and self.infos[index-1]['type'] != StateFactory.INSERT_STATE):
                print('insert')
                pos = index
                sumed_dist = {}
                for char in chars:
                    sumed_dist[char] = 0
                print(len(self.infos))
                while pos < len(self.infos) and self.infos[pos]['type'] == StateFactory.INSERT_STATE:
                    print(pos)
                    for char in chars:
                        sumed_dist[char] = sumed_dist[char] + self.infos[pos]['dist'][char]
                    pos = pos + 1
                print(sumed_dist)
            if info['type'] == StateFactory.MAIN_STATE:
                print('main')



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

st_fac = StateFactory(None, None)

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
