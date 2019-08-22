def divide(line):
    length = len(line)
    t = []
    for idx, char in enumerate(line):
        if not idx % 3 and idx + 2 < length:
            t.append(char + line[idx + 1] + line[idx + 2])
    return t


def is_illegal(list):
    reserved = ['taa', 'tag', 'tga']
    for i, tri in enumerate(list):
        if i < len(list) - 1:
            if tri in reserved:
                return True
    return False


def get_triplets(line):
    triplets = divide(line)
    while is_illegal(triplets):
        line = line[1:]
        triplets = divide(line)
    return triplets, line


def get_corrected_lines(filename):
    lines = []
    with open(filename) as file_handle:
        for i, line in enumerate(file_handle):
            curated = line.replace(',', '').replace('p', '').replace('[', '').replace(']', '').replace('\n', '')
            triplets, line = get_triplets(curated)
            lines.append(line)
    return lines


class ProbCounter:
    def __init__(self, invalid=[]):
        self.counts = {}
        self.invalid = invalid
        self.total_states = 0
        self.p = None
        self.trans = {}
        self.total_trans = 0
        self.trans_count = {
            'default': 0
        }
        self.trans_probs = None

    def add(self, key):
        self.total_trans += 1

        if key in self.trans:
            self.trans_count[self.trans[key]] += 1
        else:
            self.trans_count['default'] += 1

        if key in self.invalid:
            return
        self.total_states += 1
        if key not in self.counts:
            self.counts[key] = 1
        else:
            self.counts[key] += 1

    def calc(self):
        self.p = {}
        self.trans_probs = {}
        for key in self.trans_count:
            self.trans_probs[key] = self.trans_count[key] / self.total_trans

        for key in self.counts:
            self.p[key] = self.counts[key] / self.total_states

    def add_trans(self, key, values):
        self.trans_count[key] = 0
        for val in values:
            self.trans[val] = key


def calculate_proba(lines, order=2):
    invalid = ['a|ta', 'a|tg', 'g|ta']
    c0 = ProbCounter()
    c1 = ProbCounter()
    c2 = ProbCounter(invalid)
    c2.add_trans('end', invalid)

    for line in lines:
        for i, nu in enumerate(line):
            if i > 1:
                state_name = line[i] + '|' + line[i - 2] + line[i - 1]
                if not i % 3:
                    c0.add(state_name)
                elif i % 3 == 1:
                    c1.add(state_name)
                elif i % 3 == 2:
                    c2.add(state_name)
    c0.calc()
    c1.calc()
    c2.calc()
    return c0, c1, c2


def calculate_proba2(filename):
    invalid = ['a|ta', 'a|tg', 'g|ta']
    c0 = ProbCounter()
    c1 = ProbCounter()
    c2 = ProbCounter(invalid)
    c2.add_trans('end', invalid)
    with open(filename) as file_obj:
        for line in file_obj:
            line2 = line.replace(' ', '').replace('\n', '').lower()
            if not len(line2) % 3:
                for i, base in enumerate(line2):
                    if i > 2:
                        state_name = line2[i] + '|' + line2[i -2] + line2[i - 1]
                        if i % 3 == 0:
                            c0.add(state_name)
                        elif i % 3 == 1:
                            c1.add(state_name)
                        elif i % 3 == 2:
                            c2.add(state_name)
    c0.calc()
    c1.calc()
    c2.calc()
    return c0, c1, c2


def intron_calculator(filename):
    emission = ProbCounter()
    with open(filename) as file_obj:
        for line in file_obj:
            line2 = line.replace('\n', '').lower()
            for i, base in enumerate(line2[:-16]):
                if i > 7:
                    state = base + '|' + line2[i - 2] + line2[i - 1]
                    emission.add(state)
    emission.calc()
    return emission


def utr_exon_(filename):
    emission = ProbCounter()
    with open(filename) as file_obj:
        for line in file_obj:
            exon = line.split(' ')[0].split('P')[0].lower()
            for i, base in enumerate(exon):
                if i > 1:
                    state = base + '|' + exon[i - 2] + exon[i - 1]
                    emission.add(state)
    emission.calc()
    return emission


if __name__ == '__main__':
   utr_exon_('mcutsa.txt')