# Convert Cartesian coordinates to one or more molecular graphs

Given Cartesian coordinates in the form of a `.xyz` file, the code constructs a
list of one or more molecular graphs. In cases where there are several possible
resonance forms xyz2mol returns a list of all, otherwise just a list of one.

This code is based on the work of
DOI: [10.1002/bkcs.10334](http://dx.doi.org/10.1002/bkcs.10334)

    Yeonjoon Kim and Woo Youn Kim
    "Universal Structure Conversion Method for Organic Molecules:
    From Atomic Connectivity to Three-Dimensional Geometry"
    Bull. Korean Chem. Soc.
    2015, Vol. 36, 1769-1777

At the 2020 RDKit Virtual UGM, Jan H. Jensen presented the tool with his talk
"Dealing with organometallic molecules in RDKit"
([slides](https://github.com/rdkit/UGM_2020/blob/master/Presentations/JanJensen.pdf),
[video recording](https://www.youtube.com/watch?v=HD6IpXMVKeo)).

## Setup

Depends on `rdkit`, `numpy`, and `networkx`. Easiest to setup via anaconda/conda: 

`conda install -c conda-forge xyz2mol`

Setup for a standalone enviroment is avaliable via `Makefile`. To setup and test
simply clone the project and make.

    git clone https://github.com/jensengroup/xyz2mol

and then run the following the the `xyz2mol` folder

    make
    make test

Note, it is also possible to run the code without the `networkx` dependencies,
but is slower.


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

