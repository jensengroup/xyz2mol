# Convert Cartesian coordinates to a molecular graph

Given Cartesian coordinates in the form of a `.xyz` file, the code constructs the molecular graph.
This code is based on the work of
DOI: [10.1002/bkcs.10334](http://dx.doi.org/10.1002/bkcs.10334)

    Yeonjoon Kim and Woo Youn Kim
    "Universal Structure Conversion Method for Organic Molecules:
    From Atomic Connectivity to Three-Dimensional Geometry"
    Bull. Korean Chem. Soc.
    2015, Vol. 36, 1769-1777

## Setup

Depends on `rdkit`, `numpy`, and `networkx`. Easiest to setup via anaconda/conda.

Setup for a standalone enviroment is avaliable via `Makefile`. To setup and test simply clone the project and make.

    git clone https://github.com/jensengroup/xyz2mol

and then run the following the the `xyz2mol` folder

    make
    make test

Note, it is also possible to run the code without the `networkx` dependencies, but is slower.


## Example usage

Read in xyz file and print out the SMILES, but don't incode the chirality.

    xyz2mol.py examples/chiral_stereo_test.xyz --ignore-chiral

Read in xyz file and print out the SDF format, save it in a file

    xyz2mol.py examples/chiral_stereo_test.xyz -o sdf > save_file.sdf

Read in xyz file with a charge and print out the SMILES

    xyz2mol.py examples/acetate.xyz --charge -1

## Dependencies:

    rdkit # (version 2019.9.1 or later needed for huckel option)
    networkx

