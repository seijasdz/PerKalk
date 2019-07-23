file = 'cutsa.txt'
output = 'new_tts.exa'
before = 8
after = 8

with open(file) as file_handle, open(output, 'w') as output_handle:
    for line in file_handle:
        tokens = line.split()
        for it, token in enumerate(tokens[-1:]):
            for i, char in enumerate(token):
                if char == 'P':
                    if i + after < len(token):
                        output_handle.write(token[i - before: i + after + 1].replace('P', '').lower() + '\n')
                    else:
                        to_take = i + after - len(token) + 1
                        append = tokens[it + 1][:to_take]
                        output_handle.write(token[i - before: i + after + 1].replace('P', '').lower() + append.lower() + '\n')
