from matrix_from_aln import matrix_from_fasta
import numpy
from model_maker_utils import sequence_state_factory
from model_maker_utils import classify
from model_maker_utils import add_sequence
from model_maker_utils import spacer_states_maker
from model_maker_utils import add_variable_length_sequence
from model_maker_utils import load_long_training_examples
from model_maker_utils import get_state
from converter_to import converter_to
import calculator
from pomegranate import State
from pomegranate import HiddenMarkovModel
from pomegranate import DiscreteDistribution
from matrix_from_aln import matrix_from_exa


with open('promoter_utr_model_trained.json') as base_model_file:
    promoter_model_json = base_model_file.read()

promoter_model = HiddenMarkovModel.from_json(promoter_model_json)

matrixDonor0 = numpy.array(matrix_from_exa('new_donor0.exa'))
matrixAcceptor0 = numpy.array(matrix_from_exa('new_acceptor0.exa'))

donor0_data = classify(matrixDonor0, 2)
acceptor0_data = classify(matrixAcceptor0, 2)
no_coding_dist = calculator.intron_calculator('cuts_intron.txt').p

donor_states = sequence_state_factory(donor0_data, 'donor0')
acceptor_states = sequence_state_factory(acceptor0_data, 'acceptor0')
intron_spacer_states = spacer_states_maker(10, no_coding_dist, 'intron spacer')

utr_model = HiddenMarkovModel('utr_model')

# States
exon_state = State(DiscreteDistribution(calculator.utr_exon_('mcutsa.txt').p), name='utr exon')
intron_state = State(DiscreteDistribution(no_coding_dist), name='utr intron')

utr_model.add_model(promoter_model)
utr_model.add_state(exon_state)
utr_model.add_state(intron_state)

add_sequence(utr_model, donor_states)
add_sequence(utr_model, acceptor_states)
add_sequence(utr_model, intron_spacer_states)

utr_model.add_transition(utr_model.start, get_state(promoter_model, 'back'), 1)
utr_model.add_transition(get_state(promoter_model, 'inr7'), exon_state, 1)
utr_model.add_transition(get_state(promoter_model, 'no inr7'), exon_state, 1)


utr_model.add_transition(exon_state, exon_state, 0.60)
utr_model.add_transition(exon_state, donor_states[0], 0.000001)
utr_model.add_transition(exon_state, utr_model.end,   0.399999)

utr_model.add_transition(donor_states[-1], intron_state, 1)

utr_model.add_transition(intron_state, intron_state, 0.8)
utr_model.add_transition(intron_state, intron_spacer_states[0], 0.2)

utr_model.add_transition(intron_spacer_states[-1], acceptor_states[0], 1)

utr_model.add_transition(acceptor_states[-1], exon_state, 1)

utr_model.bake()

print(utr_model.states)

with open('utr_model_base.json', 'w',  encoding='utf-8') as out:
    out.write(utr_model.to_json())


