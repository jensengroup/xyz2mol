#
# Written by Jan H. Jensen based on this paper Yeonjoon Kim and Woo Youn Kim 
# "Universal Structure Conversion Method for Organic Molecules: From Atomic Connectivity
# to Three-Dimensional Geometry" Bull. Korean Chem. Soc. 2015, Vol. 36, 1769-1777 DOI: 10.1002/bkcs.10334
#
from rdkit import Chem
from rdkit.Chem import AllChem
import itertools
from rdkit.Chem import rdmolops
from collections import defaultdict
import copy


def getUA(maxValence_list, valence_list):
    UA = []
    DU = []
    for i, (maxValence,valence) in enumerate(zip(maxValence_list, valence_list)):
        if maxValence - valence > 0:
            UA.append(i)
            DU.append(maxValence - valence)
    return UA,DU


def get_BO(UA,AC,valences):
    BO = AC.copy()
    BO_valence = list(BO.sum(axis=1))
    UA_new, DU = getUA(valences, BO_valence)
    
    while len(UA) > 1:
        UA_pairs = list(itertools.combinations(UA, 2))

        for i,j in UA_pairs:
            if BO[i,j] > 0:
                BO[i,j] += 1
                BO[j,i] += 1
        
                BO_valence = list(BO.sum(axis=1))
                UA_new, DU_new = getUA(valences, BO_valence)
                break

        if UA_new != UA:
            UA = copy.copy(UA_new)
            DU = copy.copy(DU_new)
        else:
            break
    
    return BO, UA, DU


def BO_is_OK(BO,AC,charge,DU,atomic_valence_electrons,atomicNumList):
    q = 0
    BO_valences = list(BO.sum(axis=1))
    for i,atom in enumerate(atomicNumList):
        if atom == 1:
            continue
        else:
            q += get_atomic_charge(atom,atomic_valence_electrons[atom],BO_valences[i])
    
    if (BO-AC).sum() == sum(DU)  and charge == q:
        return True
    else:
        return False


def get_atomic_charge(atom,atomic_valence_electrons,BO_valence):
    if atom == 1:
        charge = 0
    else:
        charge = atomic_valence_electrons - 8 + BO_valence
        if atom == 16 and BO_valence == 6:
            charge = 0
        if atom == 15 and BO_valence == 5:
            charge = 0
                
    return charge


def BO2mol(mol,BO_matrix, atomicNumList,atomic_valence_electrons):
# based on code written by Paolo Toscani

    l = len(BO_matrix)
    l2 = len(atomicNumList)
    BO_valences = list(BO_matrix.sum(axis=1))
    
    if (l != l2):
        raise RuntimeError('sizes of adjMat ({0:d}) and atomicNumList '
            '{1:d} differ'.format(l, l2))

# update bond order info
    rwMol = Chem.RWMol(mol)

    bondTypeDict = {
        1: Chem.BondType.SINGLE,
        2: Chem.BondType.DOUBLE,
        3: Chem.BondType.TRIPLE
    }
    
    for i in range(l):
        for j in range(i + 1, l):
            bo = int(round(BO_matrix[i, j]))
            if (bo == 0):
                continue
            bt = bondTypeDict.get(bo, Chem.BondType.SINGLE)
            rwMol.AddBond(i, j, bt)
    mol = rwMol.GetMol()

# set atomic charges
    for i,atom in enumerate(atomicNumList):
        a = mol.GetAtomWithIdx(i)
        charge = get_atomic_charge(atom,atomic_valence_electrons[atom],BO_valences[i])
        if (abs(charge) > 0):
            a.SetFormalCharge(charge)

    rdmolops.SanitizeMol(mol)

    return mol


def AC2BO(AC,atomicNumList,charge):
    atomic_valence = defaultdict(list)
    atomic_valence[1] = [1]
    atomic_valence[6] = [4]
    atomic_valence[7] = [4,3]
    atomic_valence[8] = [2,1]
    atomic_valence[9] = [1]
    atomic_valence[14] = [4]
    atomic_valence[15] = [5,4,3]
    atomic_valence[16] = [6,4,2]
    atomic_valence[17] = [1]
    atomic_valence[35] = [1]
    atomic_valence[53] = [1]
    

    atomic_valence_electrons = defaultdict(list)
    atomic_valence_electrons[1] = 1
    atomic_valence_electrons[6] = 4
    atomic_valence_electrons[7] = 5
    atomic_valence_electrons[8] = 6
    atomic_valence_electrons[9] = 7
    atomic_valence_electrons[14] = 4
    atomic_valence_electrons[15] = 5
    atomic_valence_electrons[16] = 6
    atomic_valence_electrons[17] = 7
    atomic_valence_electrons[35] = 7
    atomic_valence_electrons[53] = 7

# make a list of valences, e.g. for CO: [[4],[2,1]]
    valences_list_of_lists = []
    for atomicNum in atomicNumList:
        valences_list_of_lists.append(atomic_valence[atomicNum])

# convert [[4],[2,1]] to [[4,2],[4,1]]
    valences_list = list(itertools.product(*valences_list_of_lists))

    best_BO = AC.copy()

# implemenation of algorithm shown in Figure 2
# UA: unsaturated atoms
# DU: degree of unsaturation (u matrix in Figure)
# best_BO: Bcurr in Figure 
#
    for valences in valences_list:
        AC_valence = list(AC.sum(axis=1))
        UA,DU = getUA(valences, AC_valence)
        DU_from_AC = copy.copy(DU)
        if len(UA) == 0:
            best_BO = AC.copy()
            break
        else:
            BO,UA,DU = get_BO(UA,AC,valences)
            if BO_is_OK(BO,AC,charge,DU_from_AC,atomic_valence_electrons,atomicNumList):
                best_BO = BO.copy()
                break
            else:
                if BO.sum() > best_BO.sum():
                    best_BO = BO.copy()

    return best_BO,atomic_valence_electrons


def AC2mol(mol,AC,atomicNumList,charge):
# convert AC matrix to bond order (BO) matrix
    BO,atomic_valence_electrons = AC2BO(AC,atomicNumList,charge)

# add BO connectivity and charge info to mol object
    mol = BO2mol(mol,BO, atomicNumList,atomic_valence_electrons)
    
    return mol


def get_proto_mol(atomicNumList):
    mol = Chem.MolFromSmarts("[#"+str(atomicNumList[0])+"]")
    rwMol = Chem.RWMol(mol)
    for i in range(1,len(atomicNumList)):
        a = Chem.Atom(atomicNumList[i])
        rwMol.AddAtom(a)
    
    mol = rwMol.GetMol()

    return mol


def xyz2AC(filename):
    import numpy as np

    symbol2number = {}
    symbol2number["H"] = 1
    symbol2number["C"] = 6
    symbol2number["N"] = 7
    symbol2number["O"] = 8
    symbol2number["F"] = 9
    symbol2number["Si"] = 14
    symbol2number["P"] = 15
    symbol2number["S"] = 16
    symbol2number["Cl"] = 17
    symbol2number["Br"] = 35
    symbol2number["I"] = 53

    atomic_symbols = []
    x_coords = []
    y_coords = []
    z_coords = []

    with open(filename, "r") as file:
        for line_number,line in enumerate(file):
            if line_number == 0:
                num_atoms = int(line)
            elif line_number == 1:
                if "charge=" in line:
                    charge = int(line.split("=")[1])
                else:
                    charge = 0
            else:
                atomic_symbol, x, y, z = line.split()
                atomic_symbols.append(atomic_symbol)
                x_coords.append(float(x))
                y_coords.append(float(y))
                z_coords.append(float(z))

    atomicNumList = []
    for symbol in atomic_symbols:
        atomicNumList.append(symbol2number[symbol])

    mol = get_proto_mol(atomicNumList)

    conf = Chem.Conformer(mol.GetNumAtoms())
    for i in range(mol.GetNumAtoms()):
        conf.SetAtomPosition(i,(x_coords[i],y_coords[i],z_coords[i]))
    mol.AddConformer(conf)

    dMat = Chem.Get3DDistanceMatrix(mol)
    pt = Chem.GetPeriodicTable()

    AC = np.zeros((num_atoms,num_atoms)).astype(int)

    for i in range(num_atoms):
        a_i = mol.GetAtomWithIdx(i)
        Rcov_i = pt.GetRcovalent(a_i.GetAtomicNum())*1.25
        for j in range(i+1,num_atoms):
            a_j = mol.GetAtomWithIdx(j)
            Rcov_j = pt.GetRcovalent(a_j.GetAtomicNum())*1.25
            if dMat[i,j] <= Rcov_i + Rcov_j:
                AC[i,j] = 1
                AC[j,i] = 1

    return AC,atomicNumList,charge,mol


def xyz2mol(filename):

# Get atom connectivity (AC) matrix, list of atomic numbers, molecular charge, 
# and mol object with no connectivity information
    AC,atomicNumList,charge,mol = xyz2AC(filename)

# Convert AC to bond order matrix and add connectivity and charge info to mol object
    new_mol = AC2mol(mol,AC,atomicNumList,charge)
    
    return new_mol


filename = "ethane.xyz"
filename = "acetate.xyz"

mol = xyz2mol(filename)

print Chem.MolToSmiles(mol)


# code to test using SMILES instead of xyz file
mol = Chem.MolFromSmiles('C=C([O-])CC')
mol = Chem.MolFromSmiles('C=C([NH3+])CC')
mol = Chem.MolFromSmiles('CC(=O)[O-]')
mol = Chem.MolFromSmiles('C[N+](=O)[O-]')
mol = Chem.MolFromSmiles('CS(CC)(=O)=O')
mol = Chem.MolFromSmiles('CS([O-])(=O)=O')
#mol = Chem.MolFromSmiles('C=C(C)CC')
#mol = Chem.MolFromSmiles('CC(C)CC')
#mol = Chem.MolFromSmiles('C=C(N)CC')
#mol = Chem.MolFromSmiles('C=C(C)C=C')
#mol = Chem.MolFromSmiles('C#CC=C')
#mol = Chem.MolFromSmiles('c1ccccc1')
#mol = Chem.MolFromSmiles('c1ccccc1c1ccccc1')

rdmolops.Kekulize(mol, clearAromaticFlags = True)
charge = Chem.GetFormalCharge(mol)
mol = Chem.AddHs(mol)
atomicNumList = [a.GetAtomicNum() for a in mol.GetAtoms()]
proto_mol = get_proto_mol(atomicNumList)

AC = Chem.GetAdjacencyMatrix(mol)

newmol = AC2mol(proto_mol,AC,atomicNumList,charge)
print Chem.MolToSmiles(newmol)


