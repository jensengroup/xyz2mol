import os
from xyz2mol import read_xyz_file, xyz2mol
from rdkit import Chem

dir_string = 'mobh35'
directory = os.fsencode(dir_string)

def canon_smiles(mol, isomeric_smiles):
    smiles = Chem.MolToSmiles(mol, isomericSmiles=isomeric_smiles)
    m = Chem.MolFromSmiles(smiles)
    smiles = Chem.MolToSmiles(m, isomericSmiles=isomeric_smiles)

    return smiles

for file in os.listdir(directory):
    filename = os.fsdecode(file)
    if 'ts' in filename or 'p' in filename:
        continue
    #print(filename)
    long_filename = '/'.join([dir_string, filename])
    atoms, charge, xyz_coordinates = read_xyz_file(long_filename)
    if '+' in filename:
        charge = 1
    mols = xyz2mol(atoms, xyz_coordinates,
        charge=charge,
        use_graph=True,
        allow_charged_fragments=True,
        embed_chiral=False,
        use_huckel=True)
    
    smiles = [canon_smiles(mol,False) for mol in mols]
    print(filename, *smiles)