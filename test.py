from rdkit import Chem
from rdkit.Chem import rdmolops
from rdkit.Chem import AllChem
import xyz2mol as x2m


if __name__ == "__main__":

    filename = "ethane.xyz"
    filename = "acetate.xyz"

    charged_fragments = True

    atomicNumList,charge,xyz_coordinates = x2m.read_xyz_file(filename)
    mol = x2m.xyz2mol(atomicNumList,charge,xyz_coordinates,charged_fragments)

    print Chem.MolToSmiles(mol)

    # code to test using SMILES instead of xyz file
    smiles_list = ['C=C([O-])CC','C=C([NH3+])CC','CC(=O)[O-]','C[N+](=O)[O-]','CS(CC)(=O)=O','CS([O-])(=O)=O',
                'C=C(C)CC', 'CC(C)CC','C=C(N)CC','C=C(C)C=C','C#CC=C','c1ccccc1','c1ccccc1c1ccccc1',
                '[NH3+]CS([O-])(=O)=O','CC(NC)=O','C[NH+]=C([O-])CC[NH+]=C([O-])C','C[NH+]=CC=C([O-])C',
                "[C+](C)(C)CC[C-](C)(C)","[CH2][CH2][CH]=[CH][CH2]"]

    for smiles in smiles_list:
        mol = Chem.MolFromSmiles(smiles)

        rdmolops.Kekulize(mol, clearAromaticFlags = True)
        charge = Chem.GetFormalCharge(mol)
        mol = Chem.AddHs(mol)
        atomicNumList = [a.GetAtomicNum() for a in mol.GetAtoms()]
        proto_mol = x2m.get_proto_mol(atomicNumList)

        AC = Chem.GetAdjacencyMatrix(mol)

        newmol = x2m.AC2mol(proto_mol,AC,atomicNumList,charge,charged_fragments)
        newmol = Chem.RemoveHs(newmol)
        newmol_smiles = Chem.MolToSmiles(newmol)

        mol = Chem.RemoveHs(mol)
        canonical_smiles = Chem.MolToSmiles(mol)
        if newmol_smiles != canonical_smiles:
            print "uh,oh", smiles, newmol_smiles

