# Decomposition of the Polarizability into Charge Transfer and Local Dipole Contributions
Computes the decomposition according to the method outlined by Laidig and Bader in J. Chem. Phys. 93, 7213–7224 (1990) (DOI: 10.1063/1.459444).
Note that these scripts do not compute the origin independent group contributuions but only the overall charge transfer and local dipole polarizabilities
of the full neutral molecule, which are origin dependent.

## Contents
1. `geometries`: Geometries of TiO$_2$ clusters
2. `submit_scripts`: includes scripts for setting up ( `setup_pol.sh`) and submitting (`submit_pol.sh`) ORCA jobs for computing the polarizabilities of the TiO$_2$ clusters via finite differences of the applied electric field. Also includes a script for calculating the atomic basins and basin moments of these clusters under different applied fields (`submit_critic2.sh`)
3. `templates`: Template input files for ORCA and Critic2 used by the scripts in `submit_scripts`.
4. `get_decomposed_polarizabilities_zero_centered.py`: python script for computing the decomposed polarizability tensors.

## Dependencies
Critic2 must be installed for computing the atomic basins. The python script requires numpy.
