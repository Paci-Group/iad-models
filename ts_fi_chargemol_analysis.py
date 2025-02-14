import pymbd as mbd
import pandas as pd
import numpy as np
import os
from gould2016 import alphamodel as am
import pickle

ang_to_bohr = 1.8897259886
gould_path = "/home/bhenderson/Desktop/atomic_polarization_project/head_benchmark_plus_clusters/python_packages/gould_bucko_2016/gould2016"
manz_ref_data = "/home/bhenderson/Desktop/atomic_polarization_project/solids/reference_data/isolated_atoms_data_with_CCSD_radial_moments_modified_headers.xlsx"
df = pd.read_excel(manz_ref_data, header=0, engine='openpyxl')    
fi_model = am.AlphaModel(file=os.path.join(gould_path, "ModelMixed_Scaled.dat"))
# ref_ion_file = "/home/bhenderson/Desktop/atomic_polarization_project/FI_references/atomic_interp_radial_moments.pkl"
ref_ion_file = "/home/bhenderson/Desktop/atomic_polarization_project/FI_references/atomic_radial_moments.pkl"
with open(ref_ion_file, 'rb') as f:
    ref_ion_moments = pickle.load(f)


def get_fi_ref_alpha(z, n, omega=0.0):
    m = np.floor(n)
    m_1 = np.ceil(n)
    f = n - m
    try:
        a_m = fi_model.GetAlpha((z,m), omega=omega)
    except Exception as e:
        print(f"WARNING: The ion Z: {z}, n_e: {m} is not defined in the reference data. Using Z: {z}, n_e: {m_1} instead.")
        m = m_1  # the charge state is not defined in the reference data. Fudge the FI approach and use integer charge instead
        a_m = fi_model.GetAlpha((z,m), omega=omega)
    try:
        a_m_1 = fi_model.GetAlpha((z, m_1), omega=omega)
    except Exception as e:
        print(f"WARNING: The ion Z: {z}, n_e: {m_1} is not defined in the reference data. Using Z: {z}, n_e: {m} instead.")
        m_1 = m  # the charge state is not defined in the reference data. Fudge the FI approach and use integer charge instead
        a_m_1 = fi_model.GetAlpha((z,m_1), omega=omega)
    return f * a_m_1 + (1 - f) * a_m

def get_fi_ref_r3(z, n):
    m = np.floor(n)
    m_1 = np.ceil(n)
    f = n - m
    return f * ref_ion_moments[(z, z, m_1)][3] + (1 - f) * ref_ion_moments[(z, z, m)][3]


def get_Rvdw(atnum, volume_ratio):
    _, _, R_vdw_ref = get_TS_ref_values(atnum) 
    return R_vdw_ref * volume_ratio**(1/3)


def get_TS_ref_values(atnum):
    alpha_0, C6, R_vdw = [mbd.vdw_params[atnum][param] for param in 'alpha_0(TS) C6(TS) R_vdw(TS)'.split()]
    return alpha_0, C6, R_vdw

def get_lattice_from_gpaw_output(filename):
    cell = []
    with open(filename, 'r') as f:
        while line := f.readline():
            if "Unit cell:" in line:
                f.readline()
                for _ in range(3):
                    line = f.readline()
                    cell.append([float(val) for val in line.split()[3:6]])
                break
    return np.array(cell)


def get_ref_r3(atnum):
    return (df[df['atomic number'] == atnum]['r^3 moment']).values[0]


def get_ref_alpha(atnum):
    return (df[df['atomic number'] == atnum]['our α']).values[0]


def get_alpha_to_r3_ratio(atnum):
    return get_ref_alpha(atnum) / get_ref_r3(atnum)


def get_aim_alpha0(atnum, r3aim, fi=False, charge=None):
    if fi:
        vratio = r3aim / get_fi_ref_r3(atnum, atnum-charge)
        alpha0 = vratio * get_fi_ref_alpha(atnum, atnum-charge)
    else:
        vratio = r3aim / get_ref_r3(atnum)
        alpha0 = vratio * get_ref_alpha(atnum)
    return alpha0, vratio


def get_all_aim_alpha0(atnums, r3aims, fi=False, charges=None):
    alpha0s, vratios = [], []
    for i, (atnum, r3aim) in enumerate(zip(atnums, r3aims)):
        if fi:
            alpha0, vratio = get_aim_alpha0(atnum, r3aim, fi, charges[i])
        else:
            alpha0, vratio = get_aim_alpha0(atnum, r3aim)
        alpha0s.append(alpha0)
        vratios.append(vratio)
    return alpha0s, vratios


def ts_analysis(ats, coords, lattice, volume_ratios, alpha0=None, idt="", beta=0.83):    
    # Get from the LibMBD Database, which uses Chu & Dalgarno reference
    alpha_0, c6_0, r_0 = mbd.from_volumes(ats, volume_ratios, kind='TS')
    if alpha0 is None:
        print(f"\n{idt}    Coordinates (Bohr) and AIM polarizabilities (a.u.) Calculated using TS-vdW Method:")
    else:    # Supplied alphas from some other reference
        print(f"\n{idt}    Coordinates (Bohr) and Supplied AIM polarizabilities (a.u.):")
        alpha_0 = alpha0
    for i, a in enumerate(alpha_0):
        print(f"{idt}       {ats[i]: <3} {coords[i,0]:13.8f} {coords[i,1]:13.8f} {coords[i,2]:13.8f} {a:13.8f}")

    # get molecular polarizability from TS polarizabilities
    iso_pol_ts = sum(alpha_0)
    print(f"\n{idt}    Isotropic Molecular Polarizability from TS-vdW: {iso_pol_ts:13.8f}")

    # get screened polarizabilities
    alpha_rsscs, alpha_rsscs_aniso, alpha_rsscs_mol, c6_rsscs, r_rsscs = mbd.screening_aim(coords, alpha_0, c6_0, r_0, beta, lattice=lattice, nfreq=1)
    print(f"\n{idt}    Coordinates (Bohr) and Isotropic AIM polarizabilities (a.u.) Calculated using TS+SCS Method:")
    for i, a in enumerate(alpha_rsscs):
        print(f"       {ats[i]: <3} {coords[i,0]:13.8f} {coords[i,1]:13.8f} {coords[i,2]:13.8f} {a:13.8f}")
    print(f"\n{idt}    Coordinates (Bohr) and Anisotropic AIM polarizabilities (a.u.) Calculated using TS+SCS Method:")
    for i, a in enumerate(alpha_rsscs_aniso):
        print(f"{idt}       {ats[i]: <3} {coords[i,0]:13.8f} {coords[i,1]:13.8f} {coords[i,2]:13.8f}")
        print(f"{idt}        " + "-"*45)
        for row in range(3):
            print(f"{idt}           {a[row,0]:13.8f} {a[row,1]:13.8f} {a[row,2]:13.8f}")
        print(f"{idt}          - Isotropic Polarizability (Tr(alpha) / 3): {np.trace(a) / 3:13.8f}\n")

    # Molecular Polarizability:
    print(f"\n{idt}    Molecular Polarizability Tensor from TS+SCS:")
    for row in range(3):
        print(f"{idt}           {alpha_rsscs_mol[row,0]:13.8f} {alpha_rsscs_mol[row,1]:13.8f} {alpha_rsscs_mol[row,2]:13.8f}")
    print(f"{idt}          - Isotropic Polarizability (Tr(alpha) / 3): {np.trace(alpha_rsscs_mol) / 3:13.8f}")
     
    # print(f"        Isotropic Sum Over Atoms: {alpha_rsscs_aniso.sum()/3:13.8f}")
    # get molecular polarizability from TS+SCS polarizabilities
    iso_pol_scs = sum(alpha_rsscs)
    print(f"\n{idt}    Isotropic Molecular Polarizability from TS+SCS: {iso_pol_scs:13.8f}\n\n")
    return iso_pol_ts, iso_pol_scs, alpha_rsscs_mol

if __name__ == "__main__":
    # Tkatchenko Scheffler 2012 for PBE TS-SCS
    # beta = 2.56
    # Ambrosetti 2014 for PBE0 MBD@rsSCS
    beta = 0.85
    lattice = None
    # loop through all materials in mbis.txt
    with open('chargemol.txt', 'r') as f:
        materials = f.read().strip().split('--')[1:]
    with open('chargemol_ts_fi_analysis.csv', 'w+') as results:
        results.write("Material,Static Polarizability (TS Manz),Static Polarizability (TS-SCS Manz),XX Manz,XY Manz,XZ Manz,YY Manz,YZ Manz,ZZ Manz,Static Polarizability (TS Chu),Static Polarizability (TS-SCS Chu),XX Chu,YY Chu,ZZ Chu,Static Polarizability (TS FI),Static Polarizability (TS-SCS FI),XX FI,XY FI,XZ FI,YY FI,YZ FI,ZZ FI\n")
        for material in materials:
            atoms = material.strip().split('\n')
            atnums = np.array([int(at.split()[2]) for at in atoms])
            atnames = np.array([at.split()[1] for at in atoms])
            coords = ang_to_bohr * np.array([[float(coord) for coord in at.split()[3:6]] for at in atoms])
            charges = np.array([float(at.split()[-2]) for at in atoms])
            aim_r3s = np.array([float(at.split()[-1]) for at in atoms])
            alphas_fi, v_ratios = get_all_aim_alpha0(atnums, aim_r3s, fi=True, charges=charges)
            alphas_fi, v_ratios = np.array(alphas_fi), np.array(v_ratios)
            
            # non-FI references
            alphas_manz, v_ratios_nonfi = get_all_aim_alpha0(atnums, aim_r3s)
            alphas_manz, v_ratios_nonfi = np.array(alphas_manz), np.array(v_ratios_nonfi)
             
            name = atoms[0].split()[0]
            print(f"Performing TS Analysis on Material {name}\n{'-'*40}\n")
            
            # get TS ad TS-SCS from Manz alphas
            print(f"*** Using the Manz CCSD(T) reference polarizabilities\n    {'-'*45}\n")
            iso_pol_ts_manz, iso_pol_scs_manz, alpha_rsscs_mol_manz = ts_analysis(atnames, coords, lattice, v_ratios_nonfi, alpha0=alphas_manz, beta=beta)
            
            # get TS and TS-SCS from Chu-Dalgarno Alphas but Manz reference r3s
            print(f"*** Using the Chu Dalgarno reference polarizabilities\n    {'-'*45}\n")
            iso_pol_ts_chu, iso_pol_scs_chu, alpha_rsscs_mol_chu = ts_analysis(atnames, coords, lattice, v_ratios_nonfi, beta=beta)

            # get TS ad TS-SCS from FI alphas
            print(f"*** Using the Gould (2016) Fractional Ion Minimal Chemistry Reference Polarizabilities\n    {'-'*45}\n")
            iso_pol_ts, iso_pol_scs, alpha_rsscs_mol = ts_analysis(atnames, coords, lattice, v_ratios, alpha0=alphas_fi, beta=beta)
            
            results.write(f"{name},{iso_pol_ts_manz:10.6f},{iso_pol_scs_manz:10.6f},{alpha_rsscs_mol_manz[0,0]:10.6f},{alpha_rsscs_mol_manz[0,1]:10.6f},{alpha_rsscs_mol_manz[0,2]:10.6f},{alpha_rsscs_mol_manz[1,1]:10.6f},{alpha_rsscs_mol_manz[1,2]:10.6f},{alpha_rsscs_mol_manz[2,2]:10.6f},{iso_pol_ts_chu:10.6f},{iso_pol_scs_chu:10.6f},{alpha_rsscs_mol_chu[0,0]:10.6f},{alpha_rsscs_mol_chu[1,1]:10.6f},{alpha_rsscs_mol_chu[2,2]:10.6f},{iso_pol_ts:10.6f},{iso_pol_scs:10.6f},{alpha_rsscs_mol[0,0]:10.6f},{alpha_rsscs_mol[0,1]:10.6f},{alpha_rsscs_mol[0,2]:10.6f},{alpha_rsscs_mol[1,1]:10.6f},{alpha_rsscs_mol[1,2]:10.6f},{alpha_rsscs_mol[2,2]:10.6f}\n")

