atg = 'new_tss.exa'
gen = 'mcutsb.txt'

with open(atg) as atg_file, open(gen) as cuts_file:
    cuts = []
    for ct in cuts_file:
        ctt = ct.lower()
        cuts.append(ctt)

    for line in atg_file:
        line2 = line.replace('\n', '')
        for ct in cuts:
            times = 0
            if line2 in ct:
                times += 1



