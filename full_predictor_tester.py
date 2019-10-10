from pathlib import Path
from xml.etree import ElementTree
from gene_ebi_to_string import to_string


def cut(tag, element, file_string):
    tokens = element.text.replace('\n', '').replace('\t', '').split('/')


def get_lines(tag, annotations, file_string):
    return [cut(tag, annotation, file_string) for annotation in annotations]


def test(tag, folder):
    path = Path(folder)
    sub_folders = [x for x in path.iterdir() if x.is_dir()]
    datas = []
    for folder in sub_folders:
        xml = folder / 'Header.ebi.xml'
        sequence_body = folder / 'sequence_body.ebi'
        if xml.exists() and sequence_body.exists():
            file_string = to_string(sequence_body)
            with xml.open() as file:
                patch = '<ROOT>' + file.read() + '</ROOT>'
            root = ElementTree.fromstring(patch)
            for ft in root.findall('FT'):
                annotations = ft.findall(tag)
                datas.append(
                    {
                        'annotations': annotations,
                        'file_string': file_string,
                        'tag': tag
                    }
                )
    lines = [get_lines(x['tag'], x['annotations'], x['file_string']) for x in datas]


if __name__ == '__main__':
    test('mRNA', '/run/media/jose/BE96A68C96A6452D/Asi/Data/')
