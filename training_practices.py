from pathlib import Path
from xml.etree import ElementTree
from gene_ebi_to_string import to_string
from pomegranate import HiddenMarkovModel
from pomegranate import State
from pomegranate import DiscreteDistribution
from converter_to import converter_to

hmmodel = HiddenMarkovModel()

back_state = State(DiscreteDistribution({'a': 0.25, 'c': 0.25, 'g': 0.25, 't': 0.25}), name='back')

fixed_state = State(DiscreteDistribution({'a': 0.45, 'c': 0.45, 'g': 0.05, 't': 0.05}), name='fixed')

hmmodel.add_state(back_state)
hmmodel.add_state(fixed_state)

hmmodel.add_transition(hmmodel.start, back_state, 1)
hmmodel.add_transition(back_state, back_state, 0.9)
hmmodel.add_transition(back_state, fixed_state, 0.1)
hmmodel.add_transition(fixed_state, fixed_state, 0.9)
hmmodel.add_transition(fixed_state, back_state, 0.1)

hmmodel.bake()

seq = list('acgtacgtaaaaccccaaa')


lopg, path = hmmodel.viterbi(seq)

print([x[1].name for x in path])

print(hmmodel.to_json())

to_fit1 = list('acgtacacacacacacac')
to_fit2 = list('acgtacgtacgtacgtacgtacgtacgtcgt')
to_fit3 = list('aaaaacccccaaacc')
to_fit4 = list('aaaaaccgcccaaaccacgtacgtacgtacgtactacgggggg')

lopg, path = hmmodel.viterbi(to_fit4)

print([x[1].name for x in path])


hmmodel.fit([to_fit1, to_fit2, to_fit3, to_fit4])

lopg, path = hmmodel.viterbi(seq)

print([x[1].name for x in path])
print(hmmodel.to_json())


