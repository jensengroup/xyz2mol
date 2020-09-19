from rdkit import Chem
from rdkit.Chem import AllChem
from collections import defaultdict
from rdkit.Chem import rdFMCS
from rdkit.Chem import rdDistGeom
from rdkit.Chem import rdEHTTools
from rdkit.Chem import rdMolDescriptors
from rdkit import rdBase

import numpy as np
from numpy import linalg as LA
from itertools import combinations

import xyz2mol

from parameters import atomic_valence_electrons, atom_aos, TMs

def get_number_of_doubly_occupied_mos(mh):
  charge = Chem.GetFormalCharge(mh)
  valence_e = [atomic_valence_electrons[a.GetAtomicNum()] for a in mh.GetAtoms()]
  number_e = sum(valence_e)-charge
  
  return number_e//2

def get_occupancy(mh,n_aos):
  '''
  Assumes singlet and doublet for even and odd number of electrons
  '''
  charge = Chem.GetFormalCharge(mh)
  valence_e = [atomic_valence_electrons[a.GetAtomicNum()] for a in mh.GetAtoms()]
  number_e = sum(valence_e)-charge
  occupancy = []
  for i in range(n_aos):
    if i < number_e//2:
      occupancy.append(2)
    elif i == number_e//2 and number_e%2 == 1:
      occupancy.append(1)
    else:
      occupancy.append(0)

  return np.array(occupancy)

def symmetrize(a):
  return a + a.T - np.diag(a.diagonal())

def get_Q_Ai(C,S,aos,occupancy):
  Q_Ai = []
  for i in range(np.count_nonzero(occupancy)):
    S_pop = occupancy[i]*np.dot(C[:,i],S[:,:])*C[:,i]
    Q_Ai.append([S_pop[ao].sum() for ao in aos])

  return np.array(Q_Ai).T

def rotate(Ci, Cj, cos, sin):
  Ci_rot =  Ci * cos + Cj * sin
  Cj_rot = -Ci * sin + Cj * cos

  return Ci_rot, Cj_rot

def get_aos(mh):
  atoms = [a.GetAtomicNum() for a in mh.GetAtoms()]
  aos = []
  count = 0
  for atom in atoms:
    aos.append(list(range(count,count+atom_aos[atom])))
    count += atom_aos[atom]

  return aos

def get_theta(C,S,aos,i,j):
  a = b = c = O0 = O1 = 0
  Spop_ii = np.dot(C[:,i],S[:,:])*C[:,i]
  Spop_jj = np.dot(C[:,j],S[:,:])*C[:,j]
  Spop_ij = 0.5 * np.dot(C[:,i],S[:,:])*C[:,j] + 0.5 * np.dot(C[:,j],S[:,:])*C[:,i]

  for atom_ao in aos:
    Aii = Spop_ii[atom_ao].sum()
    Ajj = Spop_jj[atom_ao].sum()
    Aij = Spop_ij[atom_ao].sum()

    Ad = (Aii - Ajj)
    Ao = 2.0 * Aij
    a += Ad * Ad
    b += Ao * Ao
    c += Ad * Ao

    O0 += Aij*Aij
    O1 += 0.25 * (Ajj - Aii) * (Ajj - Aii)

  Hd = a - b
  Ho = 2.0 * c
  theta = 0.5 * np.arctan2(Ho, Hd + np.sqrt(Hd * Hd + Ho * Ho))

  if abs(theta) < 1.0E-8: 
    if O1 < O0:
      theta = np.pi/4.

  return theta

def Pipek_Mezey(C,S,aos,nmos):
  V = np.copy(C)
  max_iter = 50
  for _ in range(max_iter):
    thetas = []
    mo_pairs = list(combinations(list(range(nmos)),2))
    #np.random.seed(2) #debug
    #np.random.shuffle(mo_pairs)
    for (i,j) in mo_pairs:
      theta = get_theta(V,S,aos,i,j)
      V[:,i], V[:,j] = rotate(V[:,i], V[:,j], np.cos(theta), np.sin(theta))
      thetas.append(theta)
    if max(thetas) > 0.01: #0.001:
      break

  return V

def get_MOs(mh):
  passed,res = rdEHTTools.RunMol(mh,keepOverlapAndHamiltonianMatrices=True)
  H = symmetrize(res.GetHamiltonian())
  S = symmetrize(res.GetOverlapMatrix())
  d, U = LA.eigh(S)
  S05 = U @ np.diag(np.sqrt(1/d)) @ U.T
  E, Cp = LA.eig(S05.T @H @S05)
  idx = E.argsort()[::]   
  E = E[idx]
  Cp = Cp[:,idx]
  C = S05 @ Cp

  return C, E, S

def get_TM_charges(m,charge):
  mh = Chem.Mol(m)
  mh.GetAtomWithIdx(0).SetFormalCharge(charge) #mol charge arbitrarily added to 1st atom    

  C, E, S = get_MOs(mh)
  n_aos = len(C)
  occupancy = get_occupancy(mh,n_aos)
  aos = get_aos(mh)
  number_doubly_occupied_mos = list(occupancy).count(2)
  V = Pipek_Mezey(C,S,aos,number_doubly_occupied_mos)
  Q_Ai = get_Q_Ai(V,S,aos,occupancy)  #Columns = occ MOs, rows = atoms 
  #print(Q_Ai)

  atoms = [a.GetAtomicNum() for a in mh.GetAtoms()]
  atoms2electrons = {i: 0 for i in range(len(atoms))}

  for occ,orb_pop in zip(occupancy, Q_Ai.T):
    atom_index_w_largest_pop = np.where(abs(orb_pop) == max(abs(orb_pop)))[0][0]
    #if atom_index_w_largest_pop == 14:
    #  print(atom_index_w_largest_pop,max(abs(orb_pop))) #debug
    if max(abs(orb_pop)) > 1.75:
      atoms2electrons[atom_index_w_largest_pop] += occ
      #print(atom_index_w_largest_pop,max(abs(orb_pop))) #debug

    #top_id = np.where(abs(orb_pop) > 0.5/occ)
    #if top_id[0][0] == 0: #debug
    #  print(occ, top_id, orb_pop[top_id]) #debug

  #print(atoms2electrons) #debug

  TM_charges = {}
  for i,atom in enumerate(atoms):
    if atom in TMs:
      TM_charges[i] = atomic_valence_electrons[atom]-atoms2electrons[i]

  return TM_charges

if __name__ == "__main__":
    #print(rdBase.rdkitVersion)
    file_name = 'r6.xyz'
    atoms, charge, coordinates = xyz2mol.read_xyz_file(file_name)
    #print(charge)
    AC, mol = xyz2mol.xyz2AC(atoms, coordinates, charge, use_huckel=True)
    TM_charges = get_TM_charges(mol,charge)
    print(TM_charges)