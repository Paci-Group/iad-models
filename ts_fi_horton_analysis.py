import h5py
import glob
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
ref_ion_file = "/home/bhenderson/Desktop/atomic_polarization_project/FI_references/atomic_radial_moments.pkl"
atom_db = "../horton_db_pc4"

atnum_dict = {"H": 1,
              "He": 2,
              "Li": 3,
              "Be": 4,
              "B": 5,
              "C": 6,
              "N": 7,
              "O": 8,
              "F": 9,
              "Ne": 10,
              "Na": 11,
              "Mg": 12,
              "Al": 13,
              "Si": 14,
              "P": 15,
              "S": 16,
              "Cl": 17,
              "Ar": 18,
              "K": 19,
              "Ca": 20
              }


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

def get_ref_r3(atnum):
    return (df[df['atomic number'] == atnum]['r^3 moment']).values[0]

def get_ref_moment(atom_db, atsym, part):
    at_dir = glob.glob(os.path.join(atom_db, f'*_{atsym.lower()}_*q+00', 'mult*'))[0]
    f = h5py.File(os.path.join(at_dir, f"{part}_atom_moments.h5")) 
    return f[part]['radial_moments'][0, 3]

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
    with open('horton_ts_fi_analysis.csv', 'w+') as results:
        results.write("Material,H Polarizability (TS Manz),H Polarizability (TS-SCS Manz),XX H,XY H,XZ H,YY H,YZ H,ZZ H,H Polarizability (TS FI),H Polarizability (TS-SCS FI),XX FI H,XY FI H,XZ FI H,YY FI H,YZ FI H,ZZ FI H,HI Polarizability (TS Manz),HI Polarizability (TS-SCS Manz),XX HI,XY HI,XZ HI,YY HI,YZ HI,ZZ HI,HI Polarizability (TS FI),HI Polarizability (TS-SCS FI),XX FI HI,XY FI HI,XZ FI HI,YY FI HI,YZ FI HI,ZZ FI HI,MBIS Polarizability (TS Manz),MBIS Polarizability (TS-SCS Manz),XX MBIS,XY MBIS,XZ MBIS,YY MBIS,YZ MBIS,ZZ MBIS,MBIS Polarizability (TS FI),MBIS Polarizability (TS-SCS FI),XX FI MBIS,XY FI MBIS,XZ FI MBIS,YY FI MBIS,YZ FI MBIS,ZZ FI MBIS\n")
        for filename in glob.glob("geometries/*.xyz"):
            name = os.path.basename(filename).split('.')[0]
            print(f"\nExamining Geometries for {name}:")
            print("----------------------------------")
            
            # get geometry
            with open(os.path.join(name, 'sp_property.txt'), 'r') as f:
                line = f.readline()
                while line:
                    # Atomic Symbols and Coordinates
                    if "Number of atoms" in line:
                        nat = int(line.split(':')[1])
                    elif "Coordinates" in line:
                        ats, atnums, coords = [], [], []
                        for _ in range(nat):
                            line = f.readline().split()
                            ats.append(line[1])
                            atnums.append(atnum_dict[line[1]])
                            coords.append([float(el) for el in line[2:]])
                        coords = np.array(coords) * ang_to_bohr

                    line = f.readline()


            # Retrieve MBIS Charges, Moments, and Volume Ratios
            f_mbis = h5py.File(os.path.join(name, 'mbis_output.h5'), 'r')  
            charges_mbis = f_mbis['mbis']['charges'][:]
            r3_moments_mbis = f_mbis['mbis']['radial_moments'][:, 3]
        
            # Retrieve Hirshfeld Charges, Moments, and Volume Ratios
            f_h = h5py.File(os.path.join(name, 'h_output.h5'), 'r') 
            charges_h = f_h['h']['charges'][:]
            r3_moments_h = f_h['h']['radial_moments'][:, 3]
         
            # Retrieve Iterative Hirshfeld Charges, Moments, and Volume Ratios
            f_hi = h5py.File(os.path.join(name, 'hi_output.h5'), 'r') 
            charges_hi = f_hi['hi']['charges'][:]
            r3_moments_hi = f_hi['hi']['radial_moments'][:, 3]
            
            # alphas and volume ratios using FI reference for all partitionings
            # H 
            alphas_fi_h, v_ratios_fi_h = get_all_aim_alpha0(atnums, r3_moments_h, fi=True, charges=charges_h)
            alphas_fi_h, v_ratios_fi_h = np.array(alphas_fi_h), np.array(v_ratios_fi_h)
            # HI 
            alphas_fi_hi, v_ratios_fi_hi = get_all_aim_alpha0(atnums, r3_moments_hi, fi=True, charges=charges_hi)
            alphas_fi_hi, v_ratios_fi_hi = np.array(alphas_fi_hi), np.array(v_ratios_fi_hi)
            # MBIS
            alphas_fi_mbis, v_ratios_fi_mbis = get_all_aim_alpha0(atnums, r3_moments_mbis, fi=True, charges=charges_mbis)
            alphas_fi_mbis, v_ratios_fi_mbis = np.array(alphas_fi_mbis), np.array(v_ratios_fi_mbis)
            
            # non-FI references, but still using Manz reference moments and polarizabilities
            # H
            alphas_h, v_ratios_h = get_all_aim_alpha0(atnums, r3_moments_h)
            alphas_h, v_ratios_h = np.array(alphas_h), np.array(v_ratios_h)
            # HI
            alphas_hi, v_ratios_hi = get_all_aim_alpha0(atnums, r3_moments_hi)
            alphas_hi, v_ratios_hi = np.array(alphas_hi), np.array(v_ratios_hi)
            # MBIS
            alphas_mbis, v_ratios_mbis = get_all_aim_alpha0(atnums, r3_moments_mbis)
            alphas_mbis, v_ratios_mbis = np.array(alphas_mbis), np.array(v_ratios_mbis)
             
            print(f"Performing TS Analysis on Material {name}\n{'-'*40}\n")
            
            # TS-SCS for each Partitioning Scheme
            # H
            # get TS and TS-SCS from Manz alphas
            print(f"*** Using the Hirshfeld AIM Basins \n    {'-'*45}\n")
            print(f"  +++ Using the Manz CCSD(T) reference polarizabilities\n    {'-'*45}\n")
            iso_pol_ts_h, iso_pol_scs_h, alpha_rsscs_h = ts_analysis(ats, coords, lattice, v_ratios_h, alpha0=alphas_h, beta=beta)
            # get TS and TS-SCS from FI alphas
            print(f"  +++ Using the Gould (2016) Fractional Ion Minimal Chemistry Reference Polarizabilities\n    {'-'*45}\n")
            iso_pol_ts_fi_h, iso_pol_scs_fi_h, alpha_rsscs_fi_h = ts_analysis(ats, coords, lattice, v_ratios_fi_h, alpha0=alphas_fi_h, beta=beta)
            
            # HI
            # get TS and TS-SCS from Manz alphas
            print(f"*** Using the Iterative Hirshfeld AIM Basins \n    {'-'*45}\n")
            print(f"  +++ Using the Manz CCSD(T) reference polarizabilities\n    {'-'*45}\n")
            iso_pol_ts_hi, iso_pol_scs_hi, alpha_rsscs_hi = ts_analysis(ats, coords, lattice, v_ratios_hi, alpha0=alphas_hi, beta=beta)
            # get TS and TS-SCS from FI alphas
            print(f"  +++ Using the Gould (2016) Fractional Ion Minimal Chemistry Reference Polarizabilities\n    {'-'*45}\n")
            iso_pol_ts_fi_hi, iso_pol_scs_fi_hi, alpha_rsscs_fi_hi = ts_analysis(ats, coords, lattice, v_ratios_fi_hi, alpha0=alphas_fi_hi, beta=beta)
            
            # MBIS
            # get TS and TS-SCS from Manz alphas
            print(f"*** Using the MBIS AIM Basins \n    {'-'*45}\n")
            print(f"  +++ Using the Manz CCSD(T) reference polarizabilities\n    {'-'*45}\n")
            iso_pol_ts_mbis, iso_pol_scs_mbis, alpha_rsscs_mbis = ts_analysis(ats, coords, lattice, v_ratios_mbis, alpha0=alphas_mbis, beta=beta)
            # get TS and TS-SCS from FI alphas
            print(f"  +++ Using the Gould (2016) Fractional Ion Minimal Chemistry Reference Polarizabilities\n    {'-'*45}\n")
            iso_pol_ts_fi_mbis, iso_pol_scs_fi_mbis, alpha_rsscs_fi_mbis = ts_analysis(ats, coords, lattice, v_ratios_fi_mbis, alpha0=alphas_fi_mbis, beta=beta)
            
            results.write(f"{name},{iso_pol_ts_h:10.6f},{iso_pol_scs_h:10.6f},{alpha_rsscs_h[0,0]:10.6f},{alpha_rsscs_h[0,1]:10.6f},{alpha_rsscs_h[0,2]:10.6f},{alpha_rsscs_h[1,1]:10.6f},{alpha_rsscs_h[1,2]:10.6f},{alpha_rsscs_h[2,2]:10.6f},{iso_pol_ts_fi_h:10.6f},{iso_pol_scs_fi_h:10.6f},{alpha_rsscs_fi_h[0,0]:10.6f},{alpha_rsscs_fi_h[0,1]:10.6f},{alpha_rsscs_fi_h[0,2]:10.6f},{alpha_rsscs_fi_h[1,1]:10.6f},{alpha_rsscs_fi_h[1,2]:10.6f},{alpha_rsscs_h[2,2]:10.6f},{iso_pol_ts_hi:10.6f},{iso_pol_scs_hi:10.6f},{alpha_rsscs_hi[0,0]:10.6f},{alpha_rsscs_hi[0,1]:10.6f},{alpha_rsscs_hi[0,2]:10.6f},{alpha_rsscs_hi[1,1]:10.6f},{alpha_rsscs_hi[1,2]:10.6f},{alpha_rsscs_hi[2,2]:10.6f},{iso_pol_ts_fi_hi:10.6f},{iso_pol_scs_fi_hi:10.6f},{alpha_rsscs_fi_hi[0,0]:10.6f},{alpha_rsscs_fi_hi[0,1]:10.6f},{alpha_rsscs_fi_hi[0,2]:10.6f},{alpha_rsscs_fi_hi[1,1]:10.6f},{alpha_rsscs_fi_hi[1,2]:10.6f},{alpha_rsscs_fi_hi[2,2]:10.6f},{iso_pol_ts_mbis:10.6f},{iso_pol_scs_mbis:10.6f},{alpha_rsscs_mbis[0,0]:10.6f},{alpha_rsscs_mbis[0,1]:10.6f},{alpha_rsscs_mbis[0,2]:10.6f},{alpha_rsscs_mbis[1,1]:10.6f},{alpha_rsscs_mbis[1,2]:10.6f},{alpha_rsscs_mbis[2,2]:10.6f},{iso_pol_ts_fi_mbis:10.6f},{iso_pol_scs_fi_mbis:10.6f},{alpha_rsscs_fi_mbis[0,0]:10.6f},{alpha_rsscs_fi_mbis[0,1]:10.6f},{alpha_rsscs_fi_mbis[0,2]:10.6f},{alpha_rsscs_fi_mbis[1,1]:10.6f},{alpha_rsscs_fi_mbis[1,2]:10.6f},{alpha_rsscs_fi_mbis[2,2]:10.6f}\n")

