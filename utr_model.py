from matrix_from_aln import matrix_from_fasta
import numpy
from model_maker_utils import sequence_state_factory
from model_maker_utils import classify
from model_maker_utils import add_sequence
from model_maker_utils import spacer_states_maker
from model_maker_utils import add_variable_length_sequence
from model_maker_utils import load_long_training_examples
from converter_to import converter_to
import calculator
from pomegranate import State
from pomegranate import HiddenMarkovModel
from pomegranate import DiscreteDistribution


with open('promoter_utr_model_trained.json') as base_model_file:
    promoter_model_json = base_model_file.read()

promoter_model = HiddenMarkovModel.from_json(promoter_model_json)

nvo_model = HiddenMarkovModel('nvo')
nvo_model.add_model(promoter_model)
nvo_model.add_transition(promoter_model.states[0], nvo_model.end, 1)
print(nvo_model.to_json())