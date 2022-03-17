import sys
import os
import shutil

from rdkit import Chem

sdf_file = sys.argv[1]
name = sys.argv[2]

suppl = Chem.SDMolSupplier(sdf_file)
smiles = []

for mol in suppl:
    smiles += [Chem.MolToSmiles(mol)]

if os.path.exists(name):
    shutil.rmtree(name)
os.mkdir(name)

with open(os.path.join(name, "smiles.smi"), "w") as f:
    for smi in smiles:
        f.write(smi+os.linesep)
