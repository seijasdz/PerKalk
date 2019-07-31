from pathlib import Path
from xml.etree import ElementTree
from gene_ebi_to_string import to_string


def divider(elements, sequence):
    def trim(element):
        return element.text.replace('\n', '').replace('\t', '')

    def tokenize(string):
        return string.split('/')


    el = map(trim, elements)
    el2 = map(tokenize, el)



def extractor(folder_path, lookfor):
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
            divider(elements, file_string)



if __name__ == '__main__':
    extractor(folder_path='/run/media/zippyttech/BE96A68C96A6452D/Asi/Data/', lookfor='mRNA')
