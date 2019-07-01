from collections import deque
import numpy
import pprint
from scipy.stats import chisquare
from pomegranate import HiddenMarkovModel
from pomegranate import DiscreteDistribution
from pomegranate import State

before = 3
after = 7
size = before + after
file = 'duplexW_EI'
filename = file + '.txt'
examples = 0
count_a = [0] * size
count_c = [0] * size
count_g = [0] * size
count_t = [0] * size

deck = deque(size*[0], size)
list_of_decks = []

with open(filename) as file_obj:
    with open(file + '.f', 'w') as out_file_obj:
        for line in file_obj:
            for ch in line:
                if ch != ',':
                    deck.append(ch)
                    if deck[-after] == 'p':
                        examples = examples + 1

                        llist = list(deck)
                        llist.remove('p')
                        list_of_decks.append(llist)
                        out_file_obj.write('> seq' + str(examples) + '\n')
                        for idx, val in enumerate(deck):
                            if val != 'p':
                                out_file_obj.write("%s" % val)
                            if val == 'a':
                                count_a[idx] = count_a[idx] + 1
                            elif val == 'c':
                                count_c[idx] = count_c[idx] + 1
                            elif val == 'g':
                                count_g[idx] = count_g[idx] + 1
                            elif val == 't':
                                count_t[idx] = count_t[idx] + 1
                        out_file_obj.write('\n')


og_matrix = numpy.array(list_of_decks, numpy.unicode)

final_percent = [{}] * size


MATCH = 0
NOT_MATCH = 1


def make_pair_tables(aligned_matrix, consensus):
    xy = {
        'a': 0,
        'c': 1,
        'g': 2,
        't': 3
    }
    pair_tables = {}
    for i, consensus_values in enumerate(consensus):
        for j in range(0, len(consensus)):
            if i != j:
                for row, value in enumerate(aligned_matrix[:, i]):
                    if (i, j) not in pair_tables:
                        table = numpy.array([[0, 0, 0, 0], [0, 0, 0, 0]])
                        pair_tables[(i, j)] = table
                    else:
                        table = pair_tables[(i, j)]
                    if value in consensus_values:
                        table[MATCH, xy[aligned_matrix[row, j]]] += 1
                    else:
                        table[NOT_MATCH, xy[aligned_matrix[row, j]]] += 1
    return pair_tables


def chi_square(matrix):
    rows = matrix.shape[0]
    cols = matrix.shape[1]
    expected_matrix = numpy.zeros(matrix.shape)
    suma_filas = numpy.sum(matrix, axis=0)
    suma_cols = numpy.sum(matrix, axis=1)

    for i in range(0, rows):
        for j in range(0, cols):
            expected_matrix[i, j] = (suma_filas[j] * suma_cols[i]) / numpy.sum(suma_filas)

    real_array = []
    real_expected = []

    for idx, el in enumerate(expected_matrix.ravel()):
        if el:
            real_expected.append(el)
            real_array.append(matrix.ravel()[idx])
    # print('real', real_array, real_expected)
    resu = chisquare(real_array, f_exp=real_expected)
    # print(resu)
    return resu.statistic


def chi_results(pair_tables, pos_translate):
    size = len(pos_translate)
    final_table = numpy.zeros((size, size))

    for key, array in pair_tables.items():
        # print('a table', key, array)

        if (pos_translate[key[0]] != pos_translate[key[1]])\
                and pos_translate[key[0]] != 1\
                and pos_translate[key[0]] != 2\
                and pos_translate[key[1]] != 1\
                and pos_translate[key[1]] != 2:
            # print('key:', pos_translate[key[0]], ',', pos_translate[key[1]])
            final_table[key[0], key[1]] = chi_square(array)
    return final_table


def sub_mat(matrix, idx, consensus_seq):
    total_rows = matrix.shape[0]
    in_consensus = []
    not_in_consensus = []

    for x in range(0, total_rows):
        line = matrix[x, :]
        if line[idx] in consensus_seq:
            in_consensus.append(line)
        else:
            not_in_consensus.append(line)

    return {'in_consensus': numpy.array(in_consensus, numpy.unicode),
            'not_in_consensus': numpy.array(not_in_consensus, numpy.unicode)
            }


def nucleo_percent_calc(aligned_sequences):
    counts = []
    percent = []
    total = len(aligned_sequences)
    for x in range(0, len(aligned_sequences[0])):
        counts.append({'a': 0, 'c': 0, 'g': 0, 't': 0})
    for sequence in aligned_sequences:
        for idx, val in enumerate(sequence):
            counts[idx][val] = counts[idx][val] + 1
    for count in counts:
        percent.append({
            'a': count['a'] / total,
            'c': count['c'] / total,
            'g': count['g'] / total,
            't': count['t'] / total,
        })
    print('\n\n')
    for per in percent:
        print(per)
    return percent

def find_tree(aligned_sequences, consensus_list, position_translate, min_sub_size, security):
    pair_tables = make_pair_tables(aligned_sequences, consensus_list)
    final_table = chi_results(pair_tables, position_translate)
    total_chi_square_row_sum = numpy.sum(final_table, axis=1)
    max_index = numpy.where(total_chi_square_row_sum == numpy.amax(total_chi_square_row_sum))
    # max_value = numpy.amax(total_chi_square_row_sum)
    subs = sub_mat(aligned_sequences, max_index[0][0], consensus_list[max_index[0][0]])
    # print(len(subs['not_in_consensus']), 'x', len(subs['in_consensus']))

    nucleo_percent_calc(aligned_sequences)

    this_node = {
        'leaf': True,
        'size': len(aligned_sequences),
        'index': position_translate[max_index[0][0]],
        'nucleotide': consensus_list[max_index[0][0]],
    }

    if security > 10:
        print('security!')
        return this_node

    if len(subs['not_in_consensus']) > min_sub_size and len(subs['in_consensus']) > min_sub_size:
        this_node['leaf'] = False
        this_node['not_subnode'] = find_tree(subs['not_in_consensus'], consensus_list, position_translate, min_sub_size, security + 1)
        this_node['con_subnode'] = find_tree(subs['in_consensus'], consensus_list, position_translate, min_sub_size, security + 1)

    if this_node['leaf']:
        this_node['nucleo_percent'] = nucleo_percent_calc(aligned_sequences)
    return this_node


def get_leaf_emissions(tree_node):
    if tree_node['leaf']:
        return [tree_node['nucleo_percent']]
    con_subnode = tree_node['con_subnode']
    not_subnode = tree_node['not_subnode']
    return get_leaf_emissions(con_subnode) + get_leaf_emissions(not_subnode)


def state_sequence_from(emissions, name):
    states = []
    for index, emission in enumerate(emissions):
        distribution = DiscreteDistribution(emission)
        state_name = name + '_' + str(index)
        print('creado estado', state_name)
        state = State(distribution, name=state_name)
        states.append(state)
    return states, [1] * (len(states) - 1)


def fork_sequence(model, root, branches, branch_probabilities):
    if not branch_probabilities:
        branch_probabilities = [1 / len(branches)] * len(branches)
    if len(branches) != len(branch_probabilities):
        raise Exception('lengths should be the same')
    # print('fork', branches)
    print(len(branches))
    for index, branch in enumerate(branches):
        model.add_transition(root[-1], branch[0], branch_probabilities[index])


def reunify_sequences(model, branches, merge, branch_probabilities):
    if not branch_probabilities:
        branch_probabilities = [1 / len(branches)] * len(branches)
    if len(branch_probabilities) != len(branch_probabilities):
        raise Exception('lengths should be the same')
    for index, branch in enumerate(branches):
        model.add_transition(branch[-1], merge, branch_probabilities[index])


def set_transition_probabilities(model, states, trans_probs):
    if not trans_probs:
        trans_probs = [1] * (len(states) - 1)
    for i, state in enumerate(states):
        if i < (len(states) - 1):
            model.add_transition(states[i], states[i + 1], trans_probs[i])


def add_sequence(model, states_sequence):
    for state in states_sequence:
        model.add_state(state)


donor_consensus = [['c', 'a'], ['a'], ['g'], ['g'], ['t'], ['a', 'g'], ['a'], ['g'], ['t']]
translation = [-3, -2, -1, 1, 2, 3, 4, 5, 6]
tree = find_tree(og_matrix, donor_consensus, translation, 13, 0)
pp = pprint.PrettyPrinter(indent=4)
pp.pprint(tree)

leaf_emissions = get_leaf_emissions(tree)
pp.pprint(leaf_emissions)


hm_model = HiddenMarkovModel()
mdd_states_sequences = []

for index, l_em in enumerate(leaf_emissions):
    sequence = state_sequence_from(l_em, 'donor_' + str(index))
    add_sequence(hm_model, sequence[0])
    set_transition_probabilities(hm_model, sequence[0], sequence[1])
    mdd_states_sequences.append(sequence[0])

background = State(DiscreteDistribution({'a': 0.25, 'c': 0.25, 'g': 0.25, 't': 0.25}), name='background')
hm_model.add_state(background)

hm_model.add_transition(hm_model.start, background, 0.9)
fork_sequence(hm_model, [hm_model.start], mdd_states_sequences, [0.025, 0.025, 0.025, 0.025])

hm_model.add_transition(background, background, 0.9)

fork_sequence(hm_model, [background], mdd_states_sequences, [0.025, 0.025, 0.025, 0.025])
reunify_sequences(hm_model, mdd_states_sequences, background, [1, 1, 1, 1])

hm_model.bake()
a = 'a'
c = 'c'
g = 'g'
t = 't'
seq = numpy.array([c,t,g,t,c,t,c,c,c,g,g,c,g,g,c,c,a,g,c,g,g,c,g,g,a,a,c,c,t,g,t,g,c,g,a,g,t,g,g,a,t,g,c,g,g,a,a,g,c,c,g,g,c,g,c,a,g,c,a,g,t,c,c,c,t,c,g,g,c,a,g,c,c,a,a,g,g,t,a,a,g,c,a,g,a,g,g,c,t,g,c,g,c,c,c,c,t,t,c,g,g,a,g,g,g,t,g,c,t,t,g,g,g,a,a,g,g,c,g,c,g,g,g,t,c,g,a,g,c,c,a,g,t,g,g,c,t,g,c,t,g,c,g,c,g,t,c,g])
hmm_predictions = hm_model.predict(seq)
#print("sequence: {}".format(' '.join(seq)))
#print("hmm pred: {}".format(' '.join(map( str, hmm_predictions))))
empar = []
for i, s in enumerate(seq):
    empar.append((seq[i], hm_model.states[hmm_predictions[i]].name))
print(len(hm_model.states), hm_model.states[19])
print(empar)

#for idx, val in enumerate(final_percent):
#    final_percent[idx] = {
#        'a': count_a[idx] / examples,
#        'c': count_c[idx] / examples,
#        'g': count_g[idx] / examples,
#        't': count_t[idx] / examples,
#    }

#for percent in final_percent:
#    my_sum = 0
#    for keyx, val in percent.items():
#        my_sum = my_sum + val
#    print(percent, my_sum)
