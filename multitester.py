from pathlib import Path
from xml.etree import ElementTree
from gene_ebi_to_string import to_string


def divider(elements, sequence, before, after):
    def trim(element):
        return element.text.replace('\n', '').replace('\t', '')

    def clean(sections):
        return map(
            lambda s: s.replace('join', '').replace('(', '').replace(')', ''),
            sections
        )

    def get_from(seq):
        return lambda parts: seq[int(parts[0][0]) - 1 - before: int(parts[-1][1]) + after].lower()

    trimmed_elements = map(trim, elements)
    tokens = map(lambda string: string.split('/'), trimmed_elements)
    firsts = map(lambda list_e: list_e[0], tokens)

    normals = filter(lambda s: 'complement' not in s, firsts)
    no_letters = filter(lambda s: 'A' not in s and 'Z' not in s and 'U' not in s, normals)
    clean_normals = clean(no_letters)
    pre_tuples = map(lambda line: line.split(','), clean_normals)
    tuples = map(lambda inner: list(map(lambda pre: pre.split('..'), inner)), pre_tuples)
    full_genes = map(get_from(sequence), tuples)




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
    print(list(divided))



if __name__ == '__main__':
    extract(folder_path='/run/media/zippyttech/BE96A68C96A6452D/Asi/DataEx/', lookfor='CDS', before=0, after=0)
