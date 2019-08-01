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

    def get_from(sequence):
        return lambda parts: sequence[int(parts[0][0]) - 1 - before: int(parts[-1][1]) + after]

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


def extractor(folder_path, lookfor, before, after):
    path = Path(folder_path)
    subfolders = [x for x in path.iterdir() if x.is_dir()]

    for folder in subfolders:
        xml = folder / 'Header.ebi.xml'
        sequence_body = folder / 'sequence_body.ebi'
        if xml.exists() and sequence_body.exists():
            file_string = to_string(sequence_body)
            with xml.open() as headers:
                fixed = '<ROOT>' + headers.read() + '</ROOT>'
            root = ElementTree.fromstring(fixed)
            ft = root.find('FT')
            elements = ft.findall(lookfor)
            print(divider(elements, file_string, before, after))



if __name__ == '__main__':
    extractor(folder_path='/run/media/jose/BE96A68C96A6452D/Asi/DataEx/', lookfor='CDS', before=5, after=5)
