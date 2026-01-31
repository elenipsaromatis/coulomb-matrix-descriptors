# coulomb-matrix-descriptors
This repository implements Coulomb matrix–based molecular descriptors derived from XYZ molecular geometries. The code focuses on constructing molecular representations and handling permutation dependence, without performing any supervised machine-learning training.
 

## Overview

This repository implements a small molecular representation module that converts XYZ molecular geometries into Coulomb matrix–based descriptors. The code parses `.xyz` files (element symbols + 3D coordinates), constructs the raw Coulomb matrix using nuclear charges and interatomic distances, and provides two commonly used representations to reduce sensitivity to atom ordering:

- **Eigenvalue descriptor**: a vector obtained from the eigenvalues of the Coulomb matrix.
- **Sorted Coulomb matrix**: the Coulomb matrix reordered by descending row norms.


## Methodology

For a molecule with atomic numbers `Z_i` and interatomic distances `R_ij`, the Coulomb matrix `C` is defined as:

$$
C_{ii} = 0.5 Z_i^{2.4}
$$

$$
C_{ij} = \frac{Z_i Z_j}{R_{ij}}, \quad i \neq j
$$

## Features

- Parsing of XYZ molecular structure files  
- Construction of raw Coulomb matrices  
- Eigenvalue-based permutation-invariant descriptors  
- Row-norm–sorted Coulomb matrix representations  
- Distance calculations between molecular descriptors with padding support
  

## Libraries Used
  
- NumPy   


## Input Files

This repository includes example molecular structure files in XYZ format, which are used as input for constructing Coulomb matrices:

- `H2O.xyz` – Water  
- `H2S.xyz` – Hydrogen sulfide  
- `pentane.xyz` – Pentane  

Atomic numbers are mapped from element symbols using this file:

- `atomic_number.py` – Utility for converting element symbols to nuclear charges


Developed in the context of the MSc Digital Chemistry program at Imperial College London.
