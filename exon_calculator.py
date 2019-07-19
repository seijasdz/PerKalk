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
    def __init__(self):
        self.counts = {}

    def add(self, key):
        if key not in self.counts:
            self.counts[key] = 1
        else:
            self.counts[key] += 1


def calculate_proba(lines, order=2):
    c0 = ProbCounter()
    c1 = ProbCounter()
    c2 = ProbCounter()

    for line in lines:
        for i, nu in enumerate(line):
            if i > 1:
                state_name = line[i] + '|' + line[i - 1] + line[i - 2]
                if not i % 3:
                    c0.add(state_name)
                elif i % 3 == 1:
                    c1.add(state_name)
                elif i % 3 == 2:
                    c2.add(state_name)
    print(c0.counts)


v_lines = get_corrected_lines('exonesW.txt')
calculate_proba(v_lines)