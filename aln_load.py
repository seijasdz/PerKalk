import numpy
from numpy import insert

from insertor import insert_delete_main_hmm

filename = 'duplexW_EI.aln'
to_matrix = []
with open(filename) as file_handle:
    for index, line in enumerate(file_handle):
        tokenize = line.split()
        if index > 0 and tokenize:
            row = list(tokenize[1].lower())
            print(row)
            to_matrix.append(row)

data_matrix = numpy.array(to_matrix, numpy.unicode_)

model = insert_delete_main_hmm(data_matrix)
print(model.states)

a = 'a'
c = 'c'
g = 'g'
t = 't'
seq = numpy.array(['a', 'g', 't', 't', 't', 'c', 'c', 't', 'g', 'g', 'c', 'c', 'g', 'c', 'c', 'a', 'g', 'a', 'a', 'g', 'g', 't', 'a', 't', 'g', 't'])
print(len(seq))
hmm_predictions = model.predict(seq, algorithm='viterbi')

#empar = []
#for i, s in enumerate(seq):
    #empar.append((seq[i], v_model.states[hmm_predictions[i]].name))

#print(empar)

for pre in hmm_predictions:
#    print(pre)
    print(model.states[pre].name)