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
from model_maker_utils import spacer_states_maker
from model_maker_utils import percentage_matrix_maker
from stop_example_divider import divider as stop_divider

def foo(l):
    yield from itertools.product(*([l] * 3))


c0, c1, c2 = calculator.calculate_proba2('cuts.txt')
matrixZE = numpy.array(matrix_from_exa('new_tss.exa'))

# matrixEZ = numpy.array(matrix_from_exa('new_tts.exa'))
taa_matrix, tga_matrix, tag_matrix = stop_divider('new_tts.exa')

matrixDonor0 = numpy.array(matrix_from_exa('new_donor0.exa'))
matrixDonor1 = numpy.array(matrix_from_exa('new_donor1.exa'))
matrixDonor2 = numpy.array(matrix_from_exa('new_donor2.exa'))
matrixAcceptor0 = numpy.array(matrix_from_exa('new_acceptor0.exa'))
matrixAcceptor1 = numpy.array(matrix_from_exa('new_acceptor1.exa'))
matrixAcceptor2 = numpy.array(matrix_from_exa('new_acceptor2.exa'))


polyASeqs = [
    ('AATAAA', 592),
    ('ATTAAA', 149),
    ('AGTAAA', 27),
    ('TATAAA', 32),
    ('CATAAA', 13),
    ('GATAAA', 13),
    ('AATATA', 17),
    ('AATACA', 12),
    ('AATAGA', 7),
    ('ACTAAA', 6),
    ('AAGAAA', 11),
    ('AATGAA', 8)
]

matrixPolyA = numpy.array(percentage_matrix_maker(polyASeqs))
poly_a_signal_data = classify(matrixPolyA, 2)
poly_a_states = sequence_state_factory(poly_a_signal_data, 'poly a zone ')

utr_exon_probs = calculator.utr_exon_3('mcuts.txt').p

exon3_state = State(DiscreteDistribution(utr_exon_probs), name='3utr exon')
post_poly_spacer = spacer_states_maker(15, utr_exon_probs, 'post_poly_spacer')

ze_states_data = classify(matrixZE, 2)
ze_states = sequence_state_factory(ze_states_data, 'start zone')

ez_states_taa_data = classify(numpy.array(taa_matrix), 2)
ez_states_taa = sequence_state_factory(ez_states_taa_data, 'stop zone taa')

ez_states_tga_data = classify(numpy.array(tga_matrix), 2)
ez_states_tga = sequence_state_factory(ez_states_tga_data, 'stop zone tga')

ez_states_tag_data = classify(numpy.array(tag_matrix), 2)
ez_states_tag = sequence_state_factory(ez_states_tag_data, 'stop zone tag')

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


coding_model = HiddenMarkovModel()

intron_distribution = calculator.intron_calculator('cuts_intron.txt')
back = State(DiscreteDistribution(calculator.intron_calculator('cuts_intron.txt').p), name='back')

fake_back = State(DiscreteDistribution(intron_distribution.p), name='back2')

in0 = State(DiscreteDistribution(intron_distribution.p), name='in0')
in1 = State(DiscreteDistribution(intron_distribution.p), name='in1')
in2 = State(DiscreteDistribution(intron_distribution.p), name='in2')

in0_spacers = spacer_states_maker(64, intron_distribution.p, 'in0 spacer')
in1_spacers = spacer_states_maker(64, intron_distribution.p, 'in1 spacer')
in2_spacers = spacer_states_maker(64, intron_distribution.p, 'in2 spacer')

coding_state0 = State(DiscreteDistribution(c0.p), 'coding state 0')
coding_state1 = State(DiscreteDistribution(c1.p), 'coding state 1')
coding_state2 = State(DiscreteDistribution(c2.p), 'coding state 2')

coding_model.add_state(back)
coding_model.add_state(fake_back)
coding_model.add_state(coding_state0)
coding_model.add_state(coding_state1)
coding_model.add_state(coding_state2)

coding_model.add_state(in0)
coding_model.add_state(in1)
coding_model.add_state(in2)

coding_model.add_state(exon3_state)
add_sequence(coding_model, poly_a_states)
add_sequence(coding_model, post_poly_spacer)

add_sequence(coding_model, in0_spacers)
add_sequence(coding_model, in1_spacers)
add_sequence(coding_model, in2_spacers)

add_sequence(coding_model, ze_states)

add_sequence(coding_model, ez_states_taa)
add_sequence(coding_model, ez_states_tga)
add_sequence(coding_model, ez_states_tag)

add_sequence(coding_model, donor0_states)
add_sequence(coding_model, donor1_states)
add_sequence(coding_model, donor2_states)

add_sequence(coding_model, acceptor0_states)
add_sequence(coding_model, acceptor1_states)
add_sequence(coding_model, acceptor2_states)

coding_model.add_transition(coding_model.start, back, 1.0)

coding_model.add_transition(back, back,         0.99)
coding_model.add_transition(back, ze_states[0], 0.01)

coding_model.add_transition(in0, in0,            0.99999999)
coding_model.add_transition(in0, in0_spacers[0], 0.00000001)

coding_model.add_transition(in1, in1,            0.99999999)
coding_model.add_transition(in1, in1_spacers[0], 0.00000001)

coding_model.add_transition(in2, in2,            0.99999999)
coding_model.add_transition(in2, in2_spacers[0], 0.00000001)

coding_model.add_transition(coding_state0, coding_state1, 1.0)
coding_model.add_transition(coding_state1, coding_state2, 1.0)

coding_model.add_transition(coding_state2, coding_state0,    0.9999998996999)
coding_model.add_transition(coding_state2, donor0_states[0], 0.0000000001251)
coding_model.add_transition(coding_state2, donor1_states[0], 0.0000000000800)
coding_model.add_transition(coding_state2, donor2_states[0], 0.0000000000950)
coding_model.add_transition(coding_state2, ez_states_taa[0], 0.0000000250000)
coding_model.add_transition(coding_state2, ez_states_tga[0], 0.0000000520000)
coding_model.add_transition(coding_state2, ez_states_tag[0], 0.0000000230000)

coding_model.add_transition(donor0_states[-1], in0, 1)
coding_model.add_transition(donor1_states[-1], in1, 1)
coding_model.add_transition(donor2_states[-1], in2, 1)


coding_model.add_transition(in0_spacers[-1], acceptor0_states[0], 1)
coding_model.add_transition(in1_spacers[-1], acceptor1_states[0], 1)
coding_model.add_transition(in2_spacers[-1], acceptor2_states[0], 1)

coding_model.add_transition(acceptor0_states[-1], coding_state0, 1.0)
coding_model.add_transition(acceptor1_states[-1], coding_state0, 1.0)
coding_model.add_transition(acceptor2_states[-1], coding_state0, 1.0)

coding_model.add_transition(ze_states[-1], coding_state0, 1.0)

coding_model.add_transition(ez_states_taa[-1], exon3_state, 1.0)
coding_model.add_transition(ez_states_tga[-1], exon3_state, 1.0)
coding_model.add_transition(ez_states_tag[-1], exon3_state, 1.0)

coding_model.add_transition(exon3_state, exon3_state,      0.9)
coding_model.add_transition(exon3_state, poly_a_states[0], 0.1)

coding_model.add_transition(poly_a_states[-1], post_poly_spacer[0], 1.0)
coding_model.add_transition(post_poly_spacer[-1], back, 1.0)

coding_model.bake()

with open('coding_model_base_poly.json', 'w',  encoding='utf-8') as out:
    out.write(coding_model.to_json())

