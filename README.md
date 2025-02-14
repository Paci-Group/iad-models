# Inducible Atomic Dipole Models
Python scripts for modeling molecular and crystalline geometries as an array of inducible atomic dipoles.

## Scripts
1. `ts_fi_chargemol_analysis.py`: Iterates through the molecules whose basin charges and moments are contained in the file `chargemol.txt` and computes their polarizability tensors using the volume scaling with self-consistend screening approach of Tkatchenko and Scheffler as implemented in LibMBD (more precisely the MBD@rsSCS method).
2. `ts_fi_horton_analysis.py`: Does the same thing as the above, but gets the charges and moments from `.h5` files created by HORTON. The files should each be in a directory named as the molecules they correspond with. For example, the mbis charges for BH$_3$ would be in the file `BH3/mbis_output.h5`
3. `ts_optimize_ref_plus_beta.py`: Optimizes the reference free atom polarizbilities in the model to minimize the Mean Relative Unsigned Error between the inducible atomic dipole model polarizabilities and the reference polarizabilities.
4. `decomposition/get_decomposed_polarizabilities_zero_centered.py`: Uses atomic basins from the Critic2 program to decompose the polarizability tensor into charge transfer and local dipole terms.

## Additional Data
1. `chargemol.txt`: contains the DDEC6 atomic charges and moments for the molecules from the polarizability dataset described by Hait and Head-Gordon in "How accurate are static polarizability predictions from density functional theory? An assessment over 132 species at equilibrium geometry." (2018) (DOI: 10.1039/C8CP03569E).
2. `decomposition` also contains geometries for 9 TiO$_2$ clusters and input scripts for computing atomic basins for them. Head [there](decomposition/README.md) for more information.

## Dependencies
Thes first three scripts require a slightly modified version of the libMBD package. In particular, a new function 
`screening_aim` needs to be added to the file src/pymbd/pymbd.py. The function looks like:

```python
def screening_aim(coords, alpha_0, C6, R_vdw, beta, lattice=None, nfreq=15):
    r"""Screen atomic polarizabilities.

    :param array-like coords: (a.u.) atom coordinates in rows
    :param array-like alpha_0: (a.u.) atomic polarizabilities
    :param array-like C6: (a.u.) atomic :math:`C_6` coefficients
    :param array-like R_vdw: (a.u.) atomic vdW radii
    :param float beta: MBD damping parameter :math:`\beta`
    :param array-like lattice: (a.u.) lattice vectors in rows
    :param int nfreq: number of grid points for frequency quadrature

    Returns static polarizabilities (isotropic), static atomic polarizability tensors,
    static molecular polarizability tensor, :math:`C_6` coefficients, and
    :math:`R_\mathrm{vdw}` coefficients (a.u.).
    """
    nqho = len(alpha_0)
    freq, freq_w = freq_grid(nfreq)
    omega = 4 / 3 * C6 / alpha_0**2
    alpha_dyn = [alpha_0 / (1 + (u / omega) ** 2) for u in freq]
    alpha_dyn_rsscs = []
    alpha_dyn_rsscs_aniso = []
    alpha_dyn_rsscs_mol = []
    for a in alpha_dyn:
        sigma = (np.sqrt(2 / np.pi) * a / 3) ** (1 / 3)                                               
        dipmat = dipole_matrix(                                                                       
            coords, 'fermi,dip,gg', sigma=sigma, R_vdw=R_vdw, beta=beta, lattice=lattice              
        )                                                                                             
        a_nlc = np.linalg.inv(np.diag(np.repeat(1 / a, 3)) + dipmat)                                  
        
        # atomic polarizability tensors (sum 3x3 blocks for each 3 rows)                              
        a_nlc_block = np.stack(np.split(np.stack(np.split(a_nlc, nqho, axis=0), axis=0), nqho, axis=2), axis=1)
        a_atomic = a_nlc_block.sum(axis=1)
        alpha_dyn_rsscs_aniso.append(a_atomic)                                                        

        # molecular polarizability tensor (full contraction over 3x3 blocks)                          
        a_molecular = a_atomic.sum(axis=0)
        alpha_dyn_rsscs_mol.append(a_molecular)
 
        a_contr = sum(np.sum(a_nlc[i::3, i::3], 1) for i in range(3)) / 3                             
        alpha_dyn_rsscs.append(a_contr)

    alpha_dyn_rsscs = np.stack(alpha_dyn_rsscs)                                                       
    C6_rsscs = 3 / np.pi * np.sum(freq_w[:, None] * alpha_dyn_rsscs**2, 0)                            
    R_vdw_rsscs = R_vdw * (alpha_dyn_rsscs[0, :] / alpha_0) ** (1 / 3)                                
    return alpha_dyn_rsscs[0], alpha_dyn_rsscs_aniso[0], alpha_dyn_rsscs_mol[0], C6_rsscs, R_vdw_rsscs
```

Then, the import line and the line defining the `__all__` variable need to be modified in the file `src/pymbd/__init__.py` to allow screening_aim to be importable easily from other scripts. The new lines read:

```python
from .pymbd import ang, from_volumes, mbd_energy, mbd_energy_species, screening, screening_aim
__all__ = ['mbd_energy', 'mbd_energy_species', 'screening', 'ang', 'from_volumes', 'screening_aim']
```

Finally, the scripts require the polarizability data and associated python module from the article "C6 Coefficients and Dipole Polarizabilities for All Atoms and Many Ions in Rows 1–6 of the Periodic Table" by Gould and Bucko (doi: 10.1021/acs.jctc.6b00361). That data and code was provided in the Supplementary Information of that paper. It can be installed as a python package within by placing it all within a directory gould2016. Then immediately outside of that directory, create the file  `setup.py` with the contents:

```
from setuptools import setup

if __name__ == '__main__':

    setup(
        name='gould2016',
        description='Reference free atom and ion polarizability and C6 values from DOI:10.1021/acs.jctc.6b00361',
        packages=['gould2016'],
        include_package_data=True,
        install_requires=['numpy']
        )
```

Finally from the outer directory containing `setup.py`, run `pip install -e .`.
