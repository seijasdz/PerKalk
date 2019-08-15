

text = """(329616..329646,335238..335489,336926..337210,
		338109..338262)"""


text = text.replace('(', '').replace(')', '').replace('join', '').replace('\n', '')
tokens = text.split(',')
tuples = [t.split('..') for t in tokens]
sums = [ int(t[1]) - int(t[0]) for t in tuples]
print(sums)

total = 0
for s in sums:
    total += s

print(total % 3)