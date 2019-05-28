from collections import deque
import numpy

before = 3
after = 7
size = before + after
file = 'test'
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
                        print(deck)
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
print(examples)
print(count_a)
print(count_c)
print(count_g)
print(count_t)

array_line = numpy.array(list_of_decks, numpy.unicode)
print(array_line.shape)
print(array_line)

final_percent = [{}]*size

for idx, val in enumerate(final_percent):
    final_percent[idx] = {
        'a': count_a[idx] / examples,
        'c': count_c[idx] / examples,
        'g': count_g[idx] / examples,
        't': count_t[idx] / examples,
    }

for percent in final_percent:
    my_sum = 0
    for key, val in percent.items():
        my_sum = my_sum + val
    print(percent, my_sum)
