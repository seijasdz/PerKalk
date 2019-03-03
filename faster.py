count = 0
file = 'duplexW_ZE'
filename = file + '.txt'
with open(filename) as file_obj:
    with open(file + '.f', 'w') as out_file_obj:
        for line in file_obj:
            out_file_obj.write('> seq' + str(count) + '\n')
            count = count + 1
            for ch in line:
                if ch != ',':
                    out_file_obj.write("%s" % ch)
            out_file_obj.write('\n')
