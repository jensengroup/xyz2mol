##
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


global __ATOM_LIST__
__ATOM_LIST__ = [ x.strip() for x in ['h ','he', \
      'li','be','b ','c ','n ','o ','f ','ne', \
      'na','mg','al','si','p ','s ','cl','ar', \
      'k ','ca','sc','ti','v ','cr','mn','fe','co','ni','cu', \
      'zn','ga','ge','as','se','br','kr', \
      'rb','sr','y ','zr','nb','mo','tc','ru','rh','pd','ag', \
      'cd','in','sn','sb','te','i ','xe', \
      'cs','ba','la','ce','pr','nd','pm','sm','eu','gd','tb','dy', \
      'ho','er','tm','yb','lu','hf','ta','w ','re','os','ir','pt', \
      'au','hg','tl','pb','bi','po','at','rn', \
      'fr','ra','ac','th','pa','u ','np','pu'] ]

def get_atoms_min_connectivity(AC,atomicNumList):
    heavy_atom_connectivity = []
    for atomic_number, row in zip(atomicNumList,AC):
        count = 0
        for bond, atomic_number2 in zip(row,atomicNumList):
            #print atomicNumList[i]
            if atomic_number2 != 1 and bond == 1:
                count += 1
        if atomic_number == 1:
            count = 9 # no multiple bonds to H's, so ignore
        heavy_atom_connectivity.append(count)

    min_connectivity = min(heavy_atom_connectivity)
    atoms_min_connectivity = []
    for i,connectivity in enumerate(heavy_atom_connectivity):
        if connectivity == min_connectivity:
            atoms_min_connectivity.append(i)

    return atoms_min_connectivity

def get_atom(atom):
    global __ATOM_LIST__
    atom = atom.lower()
    return __ATOM_LIST__.index(atom) + 1


def getUA(maxValence_list, valence_list):
    UA = []
    DU = []
    for i, (maxValence,valence) in enumerate(zip(maxValence_list, valence_list)):
        if maxValence - valence > 0:
            UA.append(i)
            DU.append(maxValence - valence)
    return UA,DU


def get_BO(AC,UA_try,DU,valences):
    BO = AC.copy()
    #BO_valence = list(BO.sum(axis=1))
    #UA,DU = getUA(valences, BO_valence)
    UA = list(UA_try)

    while len(DU) > 1:
        UA_pairs = itertools.combinations(UA, 2)

        for i,j in UA_pairs:
            if BO[i,j] > 0:
                BO[i,j] += 1
                BO[j,i] += 1
                break
        
        BO_valence = list(BO.sum(axis=1))
        UA_new, DU_new = getUA(valences, BO_valence)

        if DU_new != DU:
            #UA = copy.copy(UA_new)
            # change the order of UA_new elements to match that in UA_try
            UA = [UA_i for UA_i in UA_try if UA_i in UA_new]
            DU = copy.copy(DU_new)
        else:
            break

    return BO


def BO_is_OK(BO,AC,charge,DU,atomic_valence_electrons,atomicNumList,charged_fragments):

    q = 0
    if charged_fragments:
        BO_valences = list(BO.sum(axis=1))
        for i,atom in enumerate(atomicNumList):
            q += get_atomic_charge(atom,atomic_valence_electrons[atom],BO_valences[i])
            if atom == 6:
                number_of_single_bonds_to_C = list(BO[i,:]).count(1)
                if number_of_single_bonds_to_C == 2 and BO_valences[i] == 2:
                    q += 1
                if number_of_single_bonds_to_C == 3 and q + 1 < charge:
                    q += 2

    if (BO-AC).sum() == sum(DU) and charge == q:
        return True
    else:
        return False


def get_atomic_charge(atom,atomic_valence_electrons,BO_valence):
    if atom == 1:
        charge = 1 - BO_valence
    elif atom == 5:
        charge = 3 - BO_valence
    elif atom == 15 and BO_valence == 5:
        charge = 0
    elif atom == 16 and BO_valence == 6:
        charge = 0
    else:
        charge = atomic_valence_electrons - 8 + BO_valence

    return charge

def clean_charges(mol):
# this is a temporary hack. The real solution is to generate several BO matrices in AC2BO and pick the one
# with the lowest number of atomic charges
#
    #print 'clean_charges',Chem.MolToSmiles(mol)
    rxn_smarts = ['[N+:1]=[*:2]-[O-:3]>>[N+0:1]-[*:2]=[O-0:3]',
                  '[N+:1]=[*:2]-[*:3]=[*:4]-[O-:5]>>[N+0:1]-[*:2]=[*:3]-[*:4]=[O-0:5]',
                  #'[#8:1]=[#6:2][#6:3]=[#6:4][#6-:5]>>[*-:1][*:2]=[*:3][*:4]=[*+0:5]',
                  '[O:1]=[c:2][c-:3]>>[*-:1][*:2][*+0:3]',
                  '[O:1]=[C:2][C-:3]>>[*-:1][*:2]=[*+0:3]']

    fragments = Chem.GetMolFrags(mol,asMols=True)

    for i,fragment in enumerate(fragments):
        for smarts in rxn_smarts:
            patt = Chem.MolFromSmarts(smarts.split(">>")[0])
            while fragment.HasSubstructMatch(patt):
                rxn = AllChem.ReactionFromSmarts(smarts)
                ps = rxn.RunReactants((fragment,))
                fragment = ps[0][0]
        if i == 0:
            mol = fragment
        else:
            mol = Chem.CombineMols(mol,fragment)

    #print Chem.MolToSmiles(mol)

    return mol


def BO2mol(mol,BO_matrix, atomicNumList,atomic_valence_electrons,mol_charge,charged_fragments):
# based on code written by Paolo Toscani

    l = len(BO_matrix)
    l2 = len(atomicNumList)
    BO_valences = list(BO_matrix.sum(axis=1))

    if (l != l2):
        raise RuntimeError('sizes of adjMat ({0:d}) and atomicNumList '
            '{1:d} differ'.format(l, l2))

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

    if charged_fragments:
        mol = set_atomic_charges(mol,atomicNumList,atomic_valence_electrons,BO_valences,BO_matrix,mol_charge)
    else:
        mol = set_atomic_radicals(mol,atomicNumList,atomic_valence_electrons,BO_valences)

    return mol

def set_atomic_charges(mol,atomicNumList,atomic_valence_electrons,BO_valences,BO_matrix,mol_charge):
    q = 0
    for i,atom in enumerate(atomicNumList):
        a = mol.GetAtomWithIdx(i)
        charge = get_atomic_charge(atom,atomic_valence_electrons[atom],BO_valences[i])
        q += charge
        if atom == 6:
            number_of_single_bonds_to_C = list(BO_matrix[i,:]).count(1)
            if number_of_single_bonds_to_C == 2 and BO_valences[i] == 2:
                    q += 1
                    charge = 0
            if number_of_single_bonds_to_C == 3 and q + 1 < mol_charge:
                    q += 2
                    charge = 1

        if (abs(charge) > 0):
            a.SetFormalCharge(int(charge))
    #rdmolops.SanitizeMol(mol)
    mol = clean_charges(mol)

    return mol


def set_atomic_radicals(mol,atomicNumList,atomic_valence_electrons,BO_valences):
# The number of radical electrons = absolute atomic charge
    for i,atom in enumerate(atomicNumList):
        a = mol.GetAtomWithIdx(i)
        charge = get_atomic_charge(atom,atomic_valence_electrons[atom],BO_valences[i])

        if (abs(charge) > 0):
            a.SetNumRadicalElectrons(abs(charge))

    return mol


def AC2BO(AC,atomicNumList,charge,charged_fragments,quick):
    # TODO
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
    atomic_valence[32] = [4]
    atomic_valence[35] = [1]
    atomic_valence[53] = [1]


    atomic_valence_electrons = {}
    atomic_valence_electrons[1] = 1
    atomic_valence_electrons[6] = 4
    atomic_valence_electrons[7] = 5
    atomic_valence_electrons[8] = 6
    atomic_valence_electrons[9] = 7
    atomic_valence_electrons[14] = 4
    atomic_valence_electrons[15] = 5
    atomic_valence_electrons[16] = 6
    atomic_valence_electrons[17] = 7
    atomic_valence_electrons[32] = 4
    atomic_valence_electrons[35] = 7
    atomic_valence_electrons[53] = 7

# make a list of valences, e.g. for CO: [[4],[2,1]]
    valences_list_of_lists = []
    for atomicNum in atomicNumList:
        valences_list_of_lists.append(atomic_valence[atomicNum])

# convert [[4],[2,1]] to [[4,2],[4,1]]
    valences_list = list(itertools.product(*valences_list_of_lists))

    best_BO = AC.copy()

    atoms_min_connectivity = get_atoms_min_connectivity(AC,atomicNumList)

# implemenation of algorithm shown in Figure 2
# UA: unsaturated atoms
# DU: degree of unsaturation (u matrix in Figure)
# best_BO: Bcurr in Figure 
#
    is_best_BO = True
    for valences in valences_list:
        AC_valence = list(AC.sum(axis=1))
        UA,DU_from_AC = getUA(valences, AC_valence)
        if len(UA) == 0 or BO_is_OK(AC,AC,charge,DU_from_AC,atomic_valence_electrons,atomicNumList,charged_fragments):
            best_BO = AC.copy()
            break
        
        UA_atoms_min_connectivity = list(set(atoms_min_connectivity) & set(UA))
        rest = list(set(UA)-set(UA_atoms_min_connectivity))
        
        if quick:
            UA_atoms_min_connectivity_perm = [UA_atoms_min_connectivity + rest]
        else:
            UA_atoms_min_connectivity_perm = list(itertools.permutations(UA_atoms_min_connectivity))

        UA_perm = [list(a) + rest for a in UA_atoms_min_connectivity_perm]

        #print len(UA_perm)
        #print UA, UA_perm
        for UA_try in UA_perm:
            BO = get_BO(AC,UA_try,DU_from_AC,valences)
            if BO_is_OK(BO,AC,charge,DU_from_AC,atomic_valence_electrons,atomicNumList,charged_fragments):
                best_BO = BO.copy()
                is_best_BO = True
                break
            elif BO.sum() > best_BO.sum():
                    best_BO = BO.copy()
        if is_best_BO:
            break

    return best_BO,atomic_valence_electrons


def AC2mol(mol,AC,atomicNumList,charge,charged_fragments,quick):
# convert AC matrix to bond order (BO) matrix
    BO,atomic_valence_electrons = AC2BO(AC,atomicNumList,charge,charged_fragments,quick)

# add BO connectivity and charge info to mol object
    mol = BO2mol(mol,BO, atomicNumList,atomic_valence_electrons,charge,charged_fragments)

    return mol


def get_proto_mol(atomicNumList):
    mol = Chem.MolFromSmarts("[#"+str(atomicNumList[0])+"]")
    rwMol = Chem.RWMol(mol)
    for i in range(1,len(atomicNumList)):
        a = Chem.Atom(atomicNumList[i])
        rwMol.AddAtom(a)
    
    mol = rwMol.GetMol()

    return mol


def get_atomicNumList(atomic_symbols):
    atomicNumList = []
    for symbol in atomic_symbols:
        atomicNumList.append(get_atom(symbol))
    return atomicNumList


def read_xyz_file(filename):

    atomic_symbols = []
    xyz_coordinates = []

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
                xyz_coordinates.append([float(x),float(y),float(z)])

    atomicNumList = get_atomicNumList(atomic_symbols)
    
    return atomicNumList,charge,xyz_coordinates

def xyz2AC(atomicNumList,xyz):
    import numpy as np
    mol = get_proto_mol(atomicNumList)

    conf = Chem.Conformer(mol.GetNumAtoms())
    for i in range(mol.GetNumAtoms()):
        conf.SetAtomPosition(i,(xyz[i][0],xyz[i][1],xyz[i][2]))
    mol.AddConformer(conf)

    dMat = Chem.Get3DDistanceMatrix(mol)
    pt = Chem.GetPeriodicTable()

    num_atoms = len(atomicNumList)
    AC = np.zeros((num_atoms,num_atoms)).astype(int)

    for i in range(num_atoms):
        a_i = mol.GetAtomWithIdx(i)
        Rcov_i = pt.GetRcovalent(a_i.GetAtomicNum())*1.30
        for j in range(i+1,num_atoms):
            a_j = mol.GetAtomWithIdx(j)
            Rcov_j = pt.GetRcovalent(a_j.GetAtomicNum())*1.30
            if dMat[i,j] <= Rcov_i + Rcov_j:
                AC[i,j] = 1
                AC[j,i] = 1

    return AC,mol

def xyz2mol(atomicNumList,charge,xyz_coordinates,charged_fragments,quick):

# Get atom connectivity (AC) matrix, list of atomic numbers, molecular charge, 
# and mol object with no connectivity information
    AC,mol = xyz2AC(atomicNumList,xyz_coordinates)

# Convert AC to bond order matrix and add connectivity and charge info to mol object
    new_mol = AC2mol(mol,AC,atomicNumList,charge,charged_fragments,quick)

    return new_mol

if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(usage='%(prog)s [options] molecule.xyz')
    parser.add_argument('structure', metavar='structure', type=str)
    parser.add_argument('-s', '--sdf', action="store_true", help="Dump sdf file")
    args = parser.parse_args()

    
    filename = args.structure
    charged_fragments = True # alternatively radicals are made
    quick = True # try only one combination of UA atoms in AC2BO

    atomicNumList, charge, xyz_coordinates = read_xyz_file(filename)

    mol = xyz2mol(atomicNumList, charge, xyz_coordinates, charged_fragments, quick)

    if args.sdf:
        filename = filename.replace(".xyz", "")
        filename += ".sdf"
        writer = Chem.SDWriter(filename)
        writer.write(mol)

    # Canonical hack
    smiles = Chem.MolToSmiles(mol)
    m = Chem.MolFromSmiles(smiles)
    smiles = Chem.MolToSmiles(m)

    print smiles
