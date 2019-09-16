
import pytest
import rmsd
import numpy as np

from rdkit import Chem
from rdkit.Chem import rdmolops
from rdkit.Chem import AllChem

import xyz2mol as x2m

__TEST_SMILES__ = [
    'C=C([O-])CC',
    'C=C([NH3+])CC',
    'CC(=O)[O-]',
    'C[N+](=O)[O-]',
    'CS(CC)(=O)=O',
    'CS([O-])(=O)=O',
    'C=C(C)CC',
    'CC(C)CC',
    'C=C(N)CC',
    'C=C(C)C=C',
    'C#CC=C',
    'c1ccccc1',
    'c1ccccc1c1ccccc1',
    '[NH3+]CS([O-])(=O)=O',
    'CC(NC)=O',
    '[O-]c1ccccc1',
    'O=C(C=C1)C=CC1=CCC([O-])=O',
    'C#CC#C',
    'Cc1ccc(cc1)C1C=CC2C(C=CC2(C#N)C#N)=CC=1'
]

__TEST_FILES__ = [
    ("examples/ethane.xyz", 0, "CC"),
    ("examples/acetate.xyz", -1, "CC(=O)[O-]"),
    ("examples/chiral_stereo_test.xyz", 0, "C/C=C/[C@@H](C)F"),
]

def get_atoms(mol):
    atoms = [a.GetAtomicNum() for a in mol.GetAtoms()]
    return atoms

def get_mol(smiles):
    mol = Chem.MolFromSmiles(smiles)
    Chem.Kekulize(mol, clearAromaticFlags=True)
    charge = Chem.GetFormalCharge(mol)
    mol = Chem.AddHs(mol)
    return mol

def generate_structure_from_smiles(smiles):

    # Generate a 3D structure from smiles

    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)

    status = AllChem.EmbedMolecule(mol)
    status = AllChem.UFFOptimizeMolecule(mol)

    conformer = mol.GetConformer()
    coordinates = conformer.GetPositions()
    coordinates = np.array(coordinates)

    atoms = get_atoms(mol)

    return atoms, coordinates

def get_mol_skeleton(mol):

    atoms = get_atoms(mol)
    adjacent_matrix = Chem.GetAdjacencyMatrix(mol)

    return proto_mol, adjacent_matrix

@pytest.mark.parametrize("smiles", __TEST_SMILES__)
def test_smiles_from_adjacent_matrix(smiles):

    charged_fragments = True
    quick = True

    # Cut apart the smiles
    mol = get_mol(smiles)
    atoms = get_atoms(mol)
    charge = Chem.GetFormalCharge(mol)
    adjacent_matrix = Chem.GetAdjacencyMatrix(mol)

    #
    mol = Chem.RemoveHs(mol)
    canonical_smiles = Chem.MolToSmiles(mol)

    # Define new molecule template from atoms
    new_mol = x2m.get_proto_mol(atoms)

    # reconstruct the molecule from adjacent matrix, atoms and total charge
    new_mol = x2m.AC2mol(new_mol, adjacent_matrix, atoms, charge, charged_fragments, quick)
    new_mol = Chem.RemoveHs(new_mol)
    new_mol_smiles = Chem.MolToSmiles(new_mol)

    assert new_mol_smiles == canonical_smiles

    return

@pytest.mark.parametrize("smiles", __TEST_SMILES__)
def test_smiles_from_coord(smiles):

    charged_fragments = True
    quick = True

    # The answer
    mol = Chem.MolFromSmiles(smiles)
    charge = Chem.GetFormalCharge(mol)
    canonical_smiles = Chem.MolToSmiles(mol, isomericSmiles=False)

    # generate forcefield coordinates
    atoms, coordinates = generate_structure_from_smiles(smiles)

    # Generate molobj from atoms, charge and coordinates
    mol = x2m.xyz2mol(atoms, charge, coordinates, charged_fragments, quick)

    # For this test, remove chira. clean and canonical
    Chem.Kekulize(mol)
    mol = Chem.RemoveHs(mol)
    Chem.RemoveStereochemistry(mol)
    smiles = Chem.MolToSmiles(mol, isomericSmiles=False)

    # Please look away. A small hack that removes the explicit hydrogens
    mol = Chem.MolFromSmiles(smiles)
    smiles = Chem.MolToSmiles(mol)

    assert smiles == canonical_smiles

    return


@pytest.mark.parametrize("filename, charge, answer", __TEST_FILES__)
def test_smiles_from_xyz_files(filename, charge, answer):

    charged_fragments = True
    quick = True

    atoms, coordinates = rmsd.get_coordinates_xyz(filename)
    atoms = [x2m.int_atom(atom) for atom in atoms]

    mol = x2m.xyz2mol(atoms, charge, coordinates, charged_fragments, quick)
    mol = Chem.RemoveHs(mol)

    smiles = Chem.MolToSmiles(mol)

    assert smiles == answer

    return


if __name__ == "__main__":

    for smiles in __TEST_SMILES__:
        test_smiles_from_adjacent_matrix(smiles)
        print(True, smiles)

    for filename, charge, answer in __TEST_FILES__:
        test_smiles_from_xyz_files(filename, charge, answer)
        print(True, answer)

    for smiles in __TEST_SMILES__:
        test_smiles_from_coord(smiles)
        print(True, smiles)

