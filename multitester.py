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
            lambda s: s.replace('join', '').replace('(', '').replace(')', ''),
            sections
        )

    def get_from(seq, before, after):
        return lambda parts: {
            'gene': seq[parts[0][0] - before: parts[-1][1] + after].lower(),
            'parts': parts,
            'before': before,
            'after': after,
        }

    def correct_cut(cut):
        edges = cut.split('..')
        return [int(edges[0]) - 1, int(edges[1])]

    trimmed_elements = map(trim, elements)
    tokens = map(lambda string: string.split('/'), trimmed_elements)
    firsts = map(lambda list_e: list_e[0], tokens)

    normals = filter(lambda s: 'complement' not in s, firsts)
    no_letters = filter(lambda s: 'A' not in s and 'Z' not in s and 'U' not in s, normals)
    clean_normals = clean(no_letters)
    pre_tuples = map(lambda line: line.split(','), clean_normals)
    tuples = map(lambda inner: list(map(correct_cut, inner)), pre_tuples)
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

        def counter(annotated, index, count):
            if index >= len(annotated):
                return count
            if annotated[index]:
                return counter(annotated, index + 1, count + 1)
            else:
                return counter(annotated, index + 1, count)

        def validator(annotated):
            print(counter(annotated, 0, 0) / len(annotated))
        return validator

    def predict(data):
        logp, path = model.viterbi(converter_to(data['gene'], 2))
        corrected_cuts = map(corrector_generator(data['parts'][0][0], data['before']), data['parts'])
        annotated = map(annotator_generator(path), corrected_cuts)
        percent = map(validator_generator(valid_states), annotated)
        list(percent)


    return predict


if __name__ == '__main__':
    with open('hmm_model_base.json') as base_model_file:
        model_json = base_model_file.read()

    hmmodel = HiddenMarkovModel.from_json(model_json)
    genes = extract(folder_path='/run/media/jose/BE96A68C96A6452D/Asi/DataEx/', lookfor='CDS', before=50, after=10)
    valid_st = ["zone", "coding"]
    predicted = map(test(hmmodel, valid_st), genes)



    #paths = map(lambda x: x[1], predicted)
    #states = map(lambda lista_de_tuplas: list(map(lambda tupla: tupla[1].name, lista_de_tuplas)), paths)
    xerxes = list(predicted)
    print(xerxes)
