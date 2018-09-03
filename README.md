# xyz2mol

Converts an xyz file to an RDKit mol object

The code is  based on this paper Yeonjoon Kim and Woo Youn Kim "Universal Structure Conversion Method for Organic Molecules: From Atomic Connectivity to Three-Dimensional Geometry" Bull. Korean Chem. Soc. 2015, Vol. 36, 1769-1777 DOI: [10.1002/bkcs.10334](http://dx.doi.org/10.1002/bkcs.10334)

If the molecule is charged then the charge (q) must be supplied in the second line of the xyz file: charge=q= 

usage: type "python xyz2mol.py xxx.xyz"

Dependencies: [networkx](https://networkx.github.io/)

You can remove this dependency by setting quick=False and uncommenting the import statement at the top of xyz2mol.py