from rdkit import Chem
from rdkit.Chem import AllChem
from collections import defaultdict
from rdkit.Chem import rdFMCS
from rdkit.Chem import rdDistGeom
from rdkit.Chem import rdMolDescriptors
from rdkit import rdBase

from parameters import TM_charge_list, TMs
import itertools

def get_TM_charges(atoms):
    TM_charges_all = {}
    for i, atom in enumerate(atoms):
        if atom in TM_charge_list:
            TM_charges_all[i] = sorted(TM_charge_list[atom])

    keys, values = zip(*TM_charges_all.items())
    permutations_dicts = [dict(zip(keys, v)) for v in itertools.product(*values)]

    return permutations_dicts