# xyz2mol
Converts an xyz file to an RDKit mol object

The code is  based on this paper Yeonjoon Kim and Woo Youn Kim "Universal Structure Conversion Method for Organic Molecules: From Atomic Connectivity to Three-Dimensional Geometry" Bull. Korean Chem. Soc. 2015, Vol. 36, 1769-1777 DOI: 10.1002/bkcs.10334

If the molecule is charged then the charge (q) must be supplied in the second line of the xyz file: charge=q= 

usage: type python "xyz2mol.py"

I still need to add code to handle B and special rules for carbon atom charges.
Also, I only use 1 arbitrary start position when constructing the BO matrix (see Fig 1 in the paper).
