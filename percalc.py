from collections import deque
import numpy
import pprint
from scipy.stats import chisquare

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


def find_tree(aligned_sequences, consensus_list, position_translate, min_sub_size, security):
    pair_tables = make_pair_tables(aligned_sequences, consensus_list)
    final_table = chi_results(pair_tables, position_translate)
    total_chi_square_row_sum = numpy.sum(final_table, axis=1)
    max_index = numpy.where(total_chi_square_row_sum == numpy.amax(total_chi_square_row_sum))
    # max_value = numpy.amax(total_chi_square_row_sum)
    subs = sub_mat(aligned_sequences, max_index[0][0], consensus_list[max_index[0][0]])
    # print(len(subs['not_in_consensus']), 'x', len(subs['in_consensus']))


    this_node = {
        'nucleotide': consensus_list[max_index[0][0]],
    }

    if security > 10:
        print('security!')
        return this_node

    if len(subs['not_in_consensus']) > min_sub_size and len(subs['in_consensus']) > min_sub_size:
        this_node['hehe'] = 'hehe'
        this_node['not_con_subnode'] = find_tree(subs['not_in_consensus'], consensus_list, position_translate, min_sub_size, security + 1)
        this_node['con_subnode'] = find_tree(subs['in_consensus'], consensus_list, position_translate, min_sub_size, security + 1)

    return this_node

donor_consensus = [['c', 'a'], ['a'], ['g'], ['g'], ['t'], ['a', 'g'], ['a'], ['g'], ['t']]
traduc = [-3, -2, -1, 1, 2, 3, 4, 5, 6]
tree = find_tree(og_matrix, donor_consensus, traduc, 13, 0)
pp = pprint.PrettyPrinter(indent=4)
pp.pprint(tree)
for idx, val in enumerate(final_percent):
    final_percent[idx] = {
        'a': count_a[idx] / examples,
        'c': count_c[idx] / examples,
        'g': count_g[idx] / examples,
        't': count_t[idx] / examples,
    }

for percent in final_percent:
    my_sum = 0
    for keyx, val in percent.items():
        my_sum = my_sum + val
    print(percent, my_sum)
