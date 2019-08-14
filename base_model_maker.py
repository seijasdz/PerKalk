import numpy
from pomegranate import State
from pomegranate import DiscreteDistribution
from pomegranate import HiddenMarkovModel
from matrix_from_aln import matrix_from_exa
import itertools
import calculator
from model_maker_utils import sequence_state_factory
from model_maker_utils import classify
from model_maker_utils import add_sequence


def foo(l):
    yield from itertools.product(*([l] * 3))


c0, c1, c2 = calculator.calculate_proba2('cuts.txt')
matrixZE = numpy.array(matrix_from_exa('new_tss.exa'))
matrixEZ = numpy.array(matrix_from_exa('new_tts.exa'))
matrixDonor0 = numpy.array(matrix_from_exa('new_donor0.exa'))
matrixDonor1 = numpy.array(matrix_from_exa('new_donor1.exa'))
matrixDonor2 = numpy.array(matrix_from_exa('new_donor2.exa'))
matrixAcceptor0 = numpy.array(matrix_from_exa('new_acceptor0.exa'))
matrixAcceptor1 = numpy.array(matrix_from_exa('new_acceptor1.exa'))
matrixAcceptor2 = numpy.array(matrix_from_exa('new_acceptor2.exa'))


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


hmmodel = HiddenMarkovModel()

intron_distribution = calculator.intron_calculator('cuts_intron.txt')
back = State(DiscreteDistribution(intron_distribution.p), name='back')

fake_back = State(DiscreteDistribution(intron_distribution.p), name='back2')

in0 = State(DiscreteDistribution(intron_distribution.p), name='in0')
in1 = State(DiscreteDistribution(intron_distribution.p), name='in1')
in2 = State(DiscreteDistribution(intron_distribution.p), name='in2')

coding_state0 = State(DiscreteDistribution(c0.p), 'coding state 0')
coding_state1 = State(DiscreteDistribution(c1.p), 'coding state 1')
coding_state2 = State(DiscreteDistribution(c2.p), 'coding state 2')

hmmodel.add_state(back)
hmmodel.add_state(fake_back)
hmmodel.add_state(coding_state0)
hmmodel.add_state(coding_state1)
hmmodel.add_state(coding_state2)

hmmodel.add_state(in0)
hmmodel.add_state(in1)
hmmodel.add_state(in2)

hmmodel.add_transition(hmmodel.start, back, 1.0)
hmmodel.add_transition(back, back, 0.9999)

hmmodel.add_transition(in0, in0, 0.999999)
hmmodel.add_transition(in1, in1, 0.999999)
hmmodel.add_transition(in2, in2, 0.999999)

hmmodel.add_transition(coding_state0, coding_state1, 1.0)
hmmodel.add_transition(coding_state1, coding_state2, 1.0)
hmmodel.add_transition(coding_state2, coding_state0, 0.999999 - 0.00000003)

add_sequence(hmmodel, ze_states)
add_sequence(hmmodel, ez_states)

add_sequence(hmmodel, donor0_states)
add_sequence(hmmodel, donor1_states)
add_sequence(hmmodel, donor2_states)

add_sequence(hmmodel, acceptor0_states)
add_sequence(hmmodel, acceptor1_states)
add_sequence(hmmodel, acceptor2_states)

hmmodel.add_transition(coding_state2, donor0_states[0], 0.00000001)
hmmodel.add_transition(coding_state2, donor1_states[0], 0.00000001)
hmmodel.add_transition(coding_state2, donor2_states[0], 0.00000001)

hmmodel.add_transition(donor0_states[-1], in0, 1)
hmmodel.add_transition(donor1_states[-1], in1, 1)
hmmodel.add_transition(donor2_states[-1], in2, 1)

hmmodel.add_transition(in0, acceptor0_states[0], 0.000001)
hmmodel.add_transition(in1, acceptor1_states[0], 0.000001)
hmmodel.add_transition(in2, acceptor2_states[0], 0.000001)

hmmodel.add_transition(acceptor0_states[-1], coding_state0, 1.0)
hmmodel.add_transition(acceptor1_states[-1], coding_state0, 1.0)
hmmodel.add_transition(acceptor2_states[-1], coding_state0, 1.0)

hmmodel.add_transition(coding_state2, ez_states[0], 0.0000001)

hmmodel.add_transition(ez_states[-1], back, 1.0)
hmmodel.add_transition(back, ze_states[0], 0.0001)
hmmodel.add_transition(ze_states[-1], coding_state0, 1.0)

hmmodel.bake()

with open('hmm_model_base.json', 'w',  encoding='utf-8') as out:
    out.write(hmmodel.to_json())


# string = gene_ebi_to_string.to_string2('sequence_body.ebi')[365000:382042]
# print(string)
# test_seq = list(string)
# two_seq = converter_to(test_seq, 2)
# print(two_seq)
# seq = numpy.array(two_seq, numpy.unicode_)
# print(len(seq))
# logp, path = model.viterbi(seq)

# print(logp)
# count = 0
# for i, pre in enumerate(path):
    # if pre[1].name != 'back':
    #   pass
    #   print(pre[1].name, seq[i - 1])
