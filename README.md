# Convert Cartesian coordinates to one or more molecular graphs

Given Cartesian coordinates in the form of a `.xyz` file, the code constructs a list of one or more molecular graphs. In cases where there are several possible resonance forms xyz2mol returns a list of all, otherwise just a list of one.

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

## Development

Setup development enviroment with

    conda create env xyz2mol_dev -f requirements.yml
    conda activate xyz2mol_dev
    python -m pip install -r requriments.txt
    python -m pip install -r requriments.dev.txt

### Testing

Run all tests with pytest

    python -m pytest -vv test.py


### Code standard and commit

To keep code standard and ensure quality code, we use `pre-commit`.
Please install `pre-commit` in your local repository

    pre-commit install

which will install a hook into your git, and run a format and code-quality check with each commit.

### Releases

To do a release we use GitHub actions on Git Tags.
Increment the release by

    git commit . -m "We did some changes" # Do some commits
    git tag v1.1.5 # Add this tag to the newest commit
    git push origin v1.1.5 # This will trigger the release build
