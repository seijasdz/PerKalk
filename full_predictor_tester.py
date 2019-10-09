from pathlib import Path
from xml.etree import ElementTree
from gene_ebi_to_string import to_string



def test(mode, folder):
    path = Path(folder)
    sub_folders = [x for x in path.iterdir() if x.is_dir()]
    for folder in sub_folders:
        xml = folder / 'Header.ebi.xml'
        sequence_body = folder / 'sequence_body.ebi'
        if xml.exists() and sequence_body.exists():
            print(folder)
            file_string = to_string(sequence_body)
            with xml.open() as file:
                patch = '<ROOT>' + file.read() + '</ROOT>'
            root = ElementTree.fromstring(patch)
            for ft in root.findall('FT'):
                elements = ft.findall(lookfor)
                divider(elements, file_string, lookfor)
