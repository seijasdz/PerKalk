from collections import deque
before = 40
after = 35
size = before + after
file = 'TSS'
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
            for idx, ch in enumerate(line):
                # masdos = idx + 2 if idx + 2 < len(line) else idx
                # masuno = idx + 2 if idx + 1 < len(line) else idx
                # if ch != ',' and not (ch == 'a' and line[masdos] == 'p') and not (ch == 'g' and line[masuno] == 'p'):
                if ch != ',':
                    deck.append(ch)
                    if deck[-after] == 'p':
                        examples = examples + 1
                        llist = list(deck)
                        # llist.remove('p')
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