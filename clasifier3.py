import numpy
from pomegranate import State
from pomegranate import DiscreteDistribution
from pomegranate import HiddenMarkovModel
from converter_to_two import converter_to
from matrix_from_aln import matrix_from_exa
from gene_sample_extractor import seqs_from
import itertools
import exon_calculator


example = [
    ['a', 'c', 't', 'a', 't', 'g', 'c', 'a', 'a'],
    ['a', 'c', 't', 'a', 't', 'g', 'c', 'a', 'a'],
    ['a', 'c', 't', 'a', 't', 'g', 'c', 'a', 'a'],
    ['a', 'c', 't', 'a', 't', 'g', 'c', 'a', 'a'],
    ['a', 'c', 't', 'a', 't', 'g', 'c', 'a', 'a'],
    ['a', 'c', 't', 'a', 't', 'g', 'c', 'a', 'a'],
    ['a', 'c', 't', 'a', 't', 'g', 'c', 'a', 'a'],
    ['a', 'c', 't', 'a', 't', 'g', 'c', 'a', 'a'],
    ['t', 'c', 'g', 'a', 't', 'g', 'c', 'a', 'a'],
    ['t', 'c', 'g', 'a', 't', 'g', 'c', 'a', 'a'],
    ['t', 'c', 'g', 'a', 't', 'g', 'c', 'a', 'a'],
    ['t', 'c', 'g', 'a', 't', 'g', 'c', 'a', 'a'],
    ['t', 'c', 'g', 'g', 't', 'g', 'c', 'a', 'a'],
    ['t', 'c', 'g', 'a', 't', 'g', 'c', 'a', 'a'],
    ['t', 'c', 'g', 'a', 't', 'g', 'c', 'a', 'a'],

]


class HighOrderState:
    def __init__(self, columns):
        self.states_distribution = self.calculate_states(columns)

    @staticmethod
    def calculate_states(columns):
        states_distribution = {}
        for index, base1 in enumerate(columns[0]):
            name = base1 + '|'
            for i, col in enumerate(columns):
                if i:
                    name += col[index]
            if name not in states_distribution:
                states_distribution[name] = 1
            else:
                states_distribution[name] += 1

        for key, value in states_distribution.items():
            states_distribution[key] /= len(columns[0])
        return states_distribution


def classify(matrix, order=1):
    state_data = []
    for index, column in enumerate(matrix.T):
        columns = [column]
        for x in range(order, 0, -1):
            col = matrix.T[index - x]
            columns.append(col)
        if index > order - 1:
            data = HighOrderState(columns)
            state_data.append(data)
    return state_data


def sequence_state_factory(states_data, name):
    states = []
    for index, data in enumerate(states_data):
        state = State(DiscreteDistribution(data.states_distribution), name=name + str(index))
        states.append(state)
    return states


def add_sequence(model, states):
    for state in states:
        model.add_state(state)

    for index, state in enumerate(states):
        if index < len(states) - 1:
            model.add_transition(state, states[index + 1], 1.0)


def foo(l):
    yield from itertools.product(*([l] * 3))


def temporal_back(order):
    states = {}
    back_prob = 1.0 / pow(4, order + 1)
    nucleotides = ['a', 'c', 'g', 't']
    t = foo(nucleotides)
    for s in t:
        state = ''.join(s)
        states[state[:1] + '|' + state[1:]] = back_prob
    return states


v_lines = exon_calculator.get_corrected_lines('exonesW.txt')
c0, c1, c2 = exon_calculator.calculate_proba(v_lines)


matrixZE = numpy.array(matrix_from_exa('new_tss.exa'))
matrixEZ = numpy.array(matrix_from_exa('new_tts.exa'))
matrixDonor0 = numpy.array(matrix_from_exa('new_donor0.exa'))
matrixDonor1 = numpy.array(matrix_from_exa('new_donor1.exa'))
matrixDonor2 = numpy.array(matrix_from_exa('new_donor2.exa'))
matrixAcceptor0 = numpy.array(matrix_from_exa('new_acceptor0.exa'))
matrixAcceptor1 = numpy.array(matrix_from_exa('new_acceptor1.exa'))
matrixAcceptor2 = numpy.array(matrix_from_exa('new_acceptor2.exa'))

# data_matrix = numpy.array(example, numpy.unicode)

ze_states_data = classify(matrixZE, 2)
ze_states = sequence_state_factory(ze_states_data, 'start zone')

ez_states_data = classify(matrixEZ, 2)
ez_states = sequence_state_factory(ez_states_data, 'stop zone')

donor0_data = classify(matrixDonor0, 2)
donor0_states = sequence_state_factory(donor0_data, 'donor0')

donor1_data = classify(matrixDonor1, 2)
donor1_states = sequence_state_factory(donor1_data, 'donor1')

donor2_data = classify(matrixDonor2, 2)
donor2_states = sequence_state_factory(donor2_data, 'donor2')

acceptor0_data = classify(matrixAcceptor0, 2)
acceptor0_states = sequence_state_factory(acceptor0_data, 'acceptor0')

acceptor1_data = classify(matrixAcceptor1, 2)
acceptor1_states = sequence_state_factory(acceptor1_data, 'acceptor1')

acceptor2_data = classify(matrixAcceptor2, 2)
acceptor2_states = sequence_state_factory(acceptor2_data, 'acceptor2')




model = HiddenMarkovModel()

back = State(DiscreteDistribution(temporal_back(2)), name='back')

in0 = State(DiscreteDistribution(temporal_back(2)), name='in0')
in1 = State(DiscreteDistribution(temporal_back(2)), name='in1')
in2 = State(DiscreteDistribution(temporal_back(2)), name='in2')


coding_state0 = State(DiscreteDistribution(c0.p), 'coding state 0')
coding_state1 = State(DiscreteDistribution(c1.p), 'coding state 1')
coding_state2 = State(DiscreteDistribution(c2.p), 'coding state 2')

model.add_state(back)
model.add_state(coding_state0)
model.add_state(coding_state1)
model.add_state(coding_state2)

model.add_state(in0)
model.add_state(in1)
model.add_state(in2)

model.add_transition(model.start, back, 1.0)
model.add_transition(back, back, 0.9999)

model.add_transition(in0, in0, 0.9)
model.add_transition(in1, in1, 0.9)
model.add_transition(in2, in2, 0.9)

model.add_transition(coding_state0, coding_state1, 1.0)
model.add_transition(coding_state1, coding_state2, 1.0)
model.add_transition(coding_state2, coding_state0, 0.99)

add_sequence(model, ze_states)
add_sequence(model, ez_states)

add_sequence(model, donor0_states)
add_sequence(model, donor1_states)
add_sequence(model, donor2_states)

add_sequence(model, acceptor0_states)
add_sequence(model, acceptor1_states)
add_sequence(model, acceptor2_states)

model.add_transition(coding_state2, donor0_states[0], 0.003)
model.add_transition(coding_state2, donor1_states[0], 0.003)
model.add_transition(coding_state2, donor2_states[0], 0.003)

model.add_transition(donor0_states[-1], in0, 1)
model.add_transition(donor1_states[-1], in1, 1)
model.add_transition(donor2_states[-1], in2, 1)

model.add_transition(in0, acceptor0_states[0], 0.1)
model.add_transition(in1, acceptor1_states[0], 0.1)
model.add_transition(in2, acceptor2_states[0], 0.1)

model.add_transition(acceptor0_states[-1], coding_state0, 1.0)
model.add_transition(acceptor1_states[-1], coding_state0, 1.0)
model.add_transition(acceptor2_states[-1], coding_state0, 1.0)

model.add_transition(coding_state2, ez_states[0], 0.001)
model.add_transition(ez_states[-1], back, 1.0)
model.add_transition(back, ze_states[0], 0.0001)
model.add_transition(ze_states[-1], coding_state0, 1.0)

model.bake()


seqq ="AAACGTCAGCATG GTGGTATCAG CCGGCCCTTT GTCCAGCGAG AAGGCAGAGA TGAACATTCT AGAAATCAAT GAGAAATTGC GCCCCCAGTT GGCAGAGAAG  AAACAGCAGT TCAGAAACCT CAAAGAGAAA TGTTTTCTAA CTCAACTGGC CGGCTTCCTG   GCCAACCGAC AGAAGAAATA CAGTAAGATC TATAGGCTCA CCGTCATGAA AGTGATGAAT    GATGTCCTGT CTTCTCTCTG AGACACTAAA TGCTCTCTCC ATCAAAAATA ATTTCATCCT    TCCTGTACTT CTAGGAAAAC AGAAATGGGT ATTTTAACAT TTTGTTAAAG TTGGAAGACA    GAGGTACCAA AGTATTTAGC AACTTTCCAT GTTTGCAATC AGGTGGGGGT GGGACTAGAG    TTAAACTGCC ATTTATTGAT TTCTGACACA GGCACAGAAT GACCTGTTTT CTCCAAGAGG    CTCAATCATG TTTTCAAGAA TCCTCTCTGT ACCATATAAG ATCCTGCAGA CAAATAACAT"
seqq = seqq.replace(' ','').lower()
print(seqq)

back_text = 'agtagtagt'
ze_text = 'gcttacttccatgggg'
exon_text1 = 'tgcacttca'
donor0_text = 'caggtaagcag'
intron0_text = 'aaggttaaggt'
acceptor0_text = 'ggttcatatttttcaggct'
ez_text = 'gcctgatggagcct'


string = seqs_from('sequence_body.ebi')[0][0:100000]
# string = back_text + ze_text + exon_text1 + donor0_text + intron0_text + acceptor0_text + 'acgttg' + ez_text
# string = seqq
print(string)
test_seq = list(string)
two_seq = converter_to(test_seq, 2)
print(two_seq)
seq = numpy.array(two_seq, numpy.unicode_)
print(len(seq))
logp, path = model.viterbi(seq)

print(logp)
count = 0
for i, pre in enumerate(path):
     if pre[1].name != 'back':
         print(pre[1].name, seq[i - 1])