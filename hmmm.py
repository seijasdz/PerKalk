from pomegranate import *

model = HiddenMarkovModel()
#print(model.start)
s1 = State(DiscreteDistribution({'a': 0.5, 'b': 0.5}), 'st')
s2 = State(DiscreteDistribution({'a': 0.4, 'b': 0.6}), 'st2')
#print(s1, s2)
#model.add_states(s1)


model.add_transition(model.start, s1, 0.5)
model.add_transition(model.start, s2, 0.5)
#model.add_transition(s1, s1, 0.9)
#model.add_transition(s1, s2, 0.1)
#model.add_transition(s2, s1, 0.1)
#model.add_transition(s2, s2, 0.9)


model.bake()
print(model.states)