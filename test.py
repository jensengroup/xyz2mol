from rdkit import Chem
from rdkit.Chem import rdmolops
from rdkit.Chem import AllChem
import xyz2mol as x2m


if __name__ == "__main__":

    filename = "ethane.xyz"
    filename = "acetate.xyz"
    filename = "chiral_stereo_test.xyz"

    charged_fragments = True
    quick = True

    atomicNumList,charge,xyz_coordinates = x2m.read_xyz_file(filename)
    mol = x2m.xyz2mol(atomicNumList,charge,xyz_coordinates,charged_fragments,quick)

    print(Chem.MolToSmiles(mol, isomericSmiles=True))

    # code to test using SMILES instead of xyz file
    smiles_list = ['C=C([O-])CC','C=C([NH3+])CC','CC(=O)[O-]','C[N+](=O)[O-]','CS(CC)(=O)=O','CS([O-])(=O)=O',
                'C=C(C)CC', 'CC(C)CC','C=C(N)CC','C=C(C)C=C','C#CC=C','c1ccccc1','c1ccccc1c1ccccc1',
                '[NH3+]CS([O-])(=O)=O','CC(NC)=O','[O-]c1ccccc1','O=C(C=C1)C=CC1=CCC([O-])=O',
                'C#CC#C','Cc1ccc(cc1)C1C=CC2C(C=CC2(C#N)C#N)=CC=1']
    #smiles_list = ['C[NH+]=C([O-])CC[NH+]=C([O-])C','C[NH+]=CC=C([O-])C',
    #            "[C+](C)(C)CC[C-](C)(C)",'O=C(C=C1)C=CC1=CCC([O-])=O',
    #            'O=C([CH-]C=CC(C([O-])=O)=O)[O-]','[O-]c1ccccc1','CNC(C(C)=[NH+][CH-]CC(O)=O)=O',"[CH2][CH2][CH]=[CH][CH2]"]
    #smiles_list = ['Cc1ccc(cc1)C1C=CC2C(C=CC2(C#N)C#N)=CC=1']
    #smiles_list = ['CC1C=CC2C(C=CC2(C)C)=CC=1']
    #smiles_list = ['CC1=CC=C(C=CC2)C2C=C1']
    #smiles_list = ['CC1=CC=C(C2=CC=CC=C2)C=C1']
    #smiles_list = ['C1(CC2=CC=CC=C2)=CC=CC=C1']
    #smiles_list = ['[O-]c1ccccc1[O-]']
    #smiles_list = ['C[N+](=O)[O-]']
    #smiles_list = ['N#CC(C#N)=CC=C1C=CC=CC(=C1)c1ccc(cc1)[N+](=O)[O-]']
    #smiles_list = ['CNC([O-])=C([NH+]=C/CC(O)=O)C']

    for smiles in smiles_list:
        #print(smiles)
        mol = Chem.MolFromSmiles(smiles)

        Chem.Kekulize(mol, clearAromaticFlags = True)
        charge = Chem.GetFormalCharge(mol)
        mol = Chem.AddHs(mol)
        atomicNumList = [a.GetAtomicNum() for a in mol.GetAtoms()]
        proto_mol = x2m.get_proto_mol(atomicNumList)

        AC = Chem.GetAdjacencyMatrix(mol)

        newmol = x2m.AC2mol(proto_mol,AC,atomicNumList,charge,charged_fragments,quick)
        #print Chem.MolToSmiles(newmol)

        newmol = Chem.RemoveHs(newmol)
        newmol_smiles = Chem.MolToSmiles(newmol)
        #print (newmol_smiles)

        mol = Chem.RemoveHs(mol)
        canonical_smiles = Chem.MolToSmiles(mol)
        if newmol_smiles != canonical_smiles:
            print("uh,oh", smiles, newmol_smiles)

