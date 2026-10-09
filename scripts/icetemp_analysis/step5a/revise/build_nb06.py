"""Builds notebooks/06_calibration_target.ipynb (Step 5a) from the cell texts in this folder.
Run: python build_nb06.py [output.ipynb]; then execute the notebook with nbconvert."""
import sys
from pathlib import Path

import nbformat as nbf

REPO = Path('/home/jabeer/projects/glogemflow_development/glogemflow_icetemp')
OUT = Path(sys.argv[1]) if len(sys.argv) > 1 else REPO / 'notebooks' / '06_calibration_target.ipynb'
HERE = Path(__file__).resolve().parent

cells = []


def md(text):
    cells.append(nbf.v4.new_markdown_cell(text.strip('\n')))


def code(text):
    cells.append(nbf.v4.new_code_cell(text.strip('\n')))


for part in ('cells_1_intro.py', 'cells_2_data.py', 'cells_3_seasonal.py', 'cells_4_identify.py',
             'cells_5_heldout.py', 'cells_6_checks.py', 'cells_7_decision.py'):
    exec((HERE / part).read_text(encoding='utf-8'), {'md': md, 'code': code})

nb = nbf.v4.new_notebook()
nb.cells = cells
nb.metadata = {'kernelspec': {'display_name': 'Python 3', 'language': 'python', 'name': 'python3'},
               'language_info': {'name': 'python', 'pygments_lexer': 'ipython3'}}
nbf.write(nb, OUT)
print(f'wrote {OUT} ({len(cells)} cells)')
