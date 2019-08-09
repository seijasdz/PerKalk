from pathlib import Path
from xml.etree import ElementTree
from gene_ebi_to_string import to_string
from pomegranate import HiddenMarkovModel
from converter_to import converter_to


def divider(elements, sequence, before, after):
    def trim(element):
        return element.text.replace('\n', '').replace('\t', '')

    def clean(sections):
        return map(
            lambda s: (s[0].replace('join', '').replace('(', '').replace(')', ''), s[1]),
            sections
        )

    def get_from(seq, before, after):
        def annotate(data):
            parts = data[0]
            return {
                'gene': seq[parts[0][0] - before: parts[-1][1] + after].lower(),
                'parts': parts,
                'before': before,
                'after': after,
                'name': data[1]
            }
        return annotate

    def correct_cut(cut):
        edges = cut.split('..')
        return [int(edges[0]) - 1, int(edges[1])]

    trimmed_elements = map(trim, elements)
    tokens = map(lambda string: string.split('/'), trimmed_elements)
    main_data = map(lambda list_e: (list_e[0], list_e[1]), tokens)

    normals = filter(lambda s: 'complement' not in s[0], main_data)
    no_letters = filter(lambda s: 'A' not in s[0] and 'Z' not in s[0] and 'U' not in s[0], normals)
    clean_normals = clean(no_letters)
    pre_tuples = map(lambda line: (line[0].split(','), line[1]), clean_normals)
    tuples = map(lambda inner: (list(map(correct_cut, inner[0])), inner[1]), pre_tuples)
    full_genes = map(get_from(sequence, before, after), tuples)

    return list(full_genes)


def extract(folder_path, lookfor, before, after):
    def get_subfolders(path):
        return [x for x in path.iterdir() if x.is_dir()]

    def xml_string(xml):
        with xml.open() as file:
            content = '<ROOT>' + file.read() + '</ROOT>'
        return content

    def get_elements(word):
        return lambda fil: {
            'elements': ElementTree.fromstring(xml_string(fil['header'])).find('FT').findall(word),
            'sequence': to_string(fil['sequence'])
        }

    def streamline(divided):
        only_data = []
        for part in divided:
            for gene_data in part:
                only_data.append(gene_data)
        return only_data

    subfolders = get_subfolders(Path(folder_path))

    files = map(
        lambda folder: {
            'header': folder / 'Header.ebi.xml',
            'sequence': folder / 'sequence_body.ebi'
        },
        subfolders
    )

    valid_files = filter(lambda f: f['header'].exists() and f['sequence'].exists(), files)
    extracted = map(get_elements(lookfor), valid_files)
    divided = map(lambda e: divider(e['elements'], e['sequence'], before, after), extracted)
    return streamline(divided)


def test(model, valid_states):
    def corrector_generator(difference, before):
        def correct(tup):
            return tup[0] - difference + before, tup[1] - difference + before
        return correct

    def annotator_generator(path):
        def annotate(cut):
            return [x[1].name for x in path[cut[0] - 1: cut[1] - 1]]
        return annotate

    def validator_generator(valid_values):

        def counter(annotated):
            count = 0
            for char in annotated:
                if any(word in char for word in valid_values):
                    count += 1
            return count

        def validator(annotated):
            return counter(annotated) / len(annotated)
        return validator

    def predict(data):
        logp, path = model.viterbi(converter_to(data['gene'], 2))
        corrected_cuts = map(corrector_generator(data['parts'][0][0], data['before']), data['parts'])
        annotated = map(annotator_generator(path), corrected_cuts)
        percent = map(validator_generator(valid_states), annotated)
        return list(percent), data['name']

    return predict


def mean(genes):
    exon_count = 0
    tsum = 0
    checked = []
    for gene in genes:
        if gene[1] not in checked:
            checked.append(gene[1])
            tsum += sum(gene[0])
            exon_count += len(gene[0])
    print(tsum / exon_count)


if __name__ == '__main__':
    with open('hmm_model_base.json') as base_model_file:
        model_json = base_model_file.read()

    hmmodel = HiddenMarkovModel.from_json(model_json)
    genes = extract(folder_path='/run/media/jose/BE96A68C96A6452D/Asi/Data/', lookfor='CDS', before=50, after=10)
    valid_st = ["zone", 'coding',
                'acceptor016', 'acceptor017', 'acceptor018',
                'acceptor116', 'acceptor117', 'acceptor118', 'acceptor119', 'acceptor120',
                'acceptor216', 'acceptor217', 'acceptor218', 'acceptor219',
                'donor00', 'donor01', 'donor02',
                'donor10', 'donor11', 'donor12', 'donor13',
                'donor20', 'donor21', 'donor22', 'donor23', 'donor24',
                ]
    predicted = map(test(hmmodel, valid_st), genes)

    exon_match = list(predicted)
    # print(exon_match)
    mean(exon_match)
