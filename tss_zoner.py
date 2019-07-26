atg = 'new_tss.exa'
gen = 'mcutsb.txt'

with open(atg) as atg_file, open(gen) as cuts_file:
    cuts = []
    for ct in cuts_file:
        ctt = ct.lower()
        cuts.append(ctt)

    for line in atg_file:
        line2 = line.replace('\n', '')
        found = False
        for ct in cuts:
            ct2 = ct.replace('p', '')
            times = 0
            if line2 in ct2:
                found = True
                print(ct.index('p'))
                print(ct2[25: ct2.index(line2) + len(line2)-6])
        if not found:
            pass
            # print(line2)


