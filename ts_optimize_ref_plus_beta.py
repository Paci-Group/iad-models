import pymbd as mbd
import pandas as pd
import numpy as np
import os
from scipy.optimize import minimize

ang_to_bohr = 1.8897259886
manz_ref_data = "/home/bhenderson/Desktop/atomic_polarization_project/solids/reference_data/isolated_atoms_data_with_CCSD_radial_moments_modified_headers.xlsx"
df = pd.read_excel(manz_ref_data, header=0, engine='openpyxl')    

def get_Rvdw(atnum, volume_ratio):
    _, _, R_vdw_ref = get_TS_ref_values(atnum) 
    return R_vdw_ref * volume_ratio**(1/3)


def get_TS_ref_values(atnum):
    alpha_0, C6, R_vdw = [mbd.vdw_params[atnum][param] for param in 'alpha_0(TS) C6(TS) R_vdw(TS)'.split()]
    return alpha_0, C6, R_vdw


def get_ref_r3(atnum):
    return (df[df['atomic number'] == atnum]['r^3 moment']).values[0]


def get_ref_alpha(atnum):
    return (df[df['atomic number'] == atnum]['our α']).values[0]


def get_alpha_to_r3_ratio(atnum):
    return get_ref_alpha(atnum) / get_ref_r3(atnum)


def get_aim_alpha0(atnum, r3aim, alpha_ref=None):
    vratio = r3aim / get_ref_r3(atnum)
    if alpha_refs is None:
        alpha0 = vratio * get_ref_alpha(atnum)
    else:
        alpha0 = vratio * alpha_ref

    return alpha0, vratio


def get_all_aim_alpha0(atnums, r3aims, alpha_refs=None, alpha_refs_key=None):
    alpha0s, vratios = [], []
    for i, (atnum, r3aim) in enumerate(zip(atnums, r3aims)):
        alpha_ref = alpha_refs[alpha_refs_key.index(atnum)]
        alpha0, vratio = get_aim_alpha0(atnum, r3aim, alpha_ref=alpha_ref)
        alpha0s.append(alpha0)
        vratios.append(vratio)
    return alpha0s, vratios


def ts_lite(ats, coords, lattice, volume_ratios, alpha_0, beta=0.83):    
    # Get from the LibMBD Database, which uses Chu & Dalgarno reference
    _, c6_0, r_0 = mbd.from_volumes(ats, volume_ratios, kind='TS')
    # get screened polarizabilities
    alpha_rsscs, alpha_rsscs_aniso, alpha_rsscs_mol, c6_rsscs, r_rsscs = mbd.screening_aim(coords, alpha_0, c6_0, r_0, beta, lattice=lattice, nfreq=1)
    # get molecular polarizability from TS+SCS polarizabilities
    iso_pol_scs = sum(alpha_rsscs)
    return iso_pol_scs, alpha_rsscs_mol


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


def objective_mre_plus_beta(alpha_refs_plus_beta, clusters, dft_isos, alpha_refs_key):
    tot_rue = 0
    alpha_refs, beta = alpha_refs_plus_beta[:-1], alpha_refs_plus_beta[-1]

    # for each cluster
    for cluster, dft_iso in zip(clusters, dft_isos):

        # get the ts alpha_iso with the given free atom references
        atnums = cluster['atnums']
        atnames = cluster['atnames']
        coords = cluster['coords']
        aim_r3s = cluster['r3s']

        # non-FI starting guess
        alphas_new, v_ratios_nonfi = get_all_aim_alpha0(atnums, aim_r3s, alpha_refs, alpha_refs_key=alpha_refs_key)
        alphas_new, v_ratios_nonfi = np.array(alphas_new), np.array(v_ratios_nonfi)


        # calculate the relative error of the alpha_iso
        iso, alpha_aniso = ts_lite(atnames, coords, None, v_ratios_nonfi, alphas_new, beta=beta) 
        rue = np.abs(iso - dft_iso) / dft_iso

        # add to total errors
        tot_rue +=rue

    # return mean relative (absolute) error
    return tot_rue / len(clusters)

def objective_mrse_plus_beta(alpha_refs_plus_beta, clusters, dft_isos, alpha_refs_key):
    tot_rue = 0
    alpha_refs, beta = alpha_refs_plus_beta[:-1], alpha_refs_plus_beta[-1]

    # for each cluster
    for cluster, dft_iso in zip(clusters, dft_isos):

        # get the ts alpha_iso with the given free atom references
        atnums = cluster['atnums']
        atnames = cluster['atnames']
        coords = cluster['coords']
        aim_r3s = cluster['r3s']

        # non-FI starting guess
        alphas_new, v_ratios_nonfi = get_all_aim_alpha0(atnums, aim_r3s, alpha_refs, alpha_refs_key=alpha_refs_key)
        alphas_new, v_ratios_nonfi = np.array(alphas_new), np.array(v_ratios_nonfi)


        # calculate the relative error of the alpha_iso
        iso, alpha_aniso = ts_lite(atnames, coords, None, v_ratios_nonfi, alphas_new, beta=beta) 
        rue = (iso - dft_iso) / dft_iso

        # add to total errors
        tot_rue +=rue

    # return mean relative (absolute) error
    return tot_rue / len(clusters)

def objective_mre(alpha_refs, clusters, dft_isos, alpha_refs_key, beta):
    tot_rue = 0

    # for each cluster
    for cluster, dft_iso in zip(clusters, dft_isos):

        # get the ts alpha_iso with the given free atom references
        atnums = cluster['atnums']
        atnames = cluster['atnames']
        coords = cluster['coords']
        aim_r3s = cluster['r3s']

        # non-FI starting guess
        alphas_new, v_ratios_nonfi = get_all_aim_alpha0(atnums, aim_r3s, alpha_refs, alpha_refs_key=alpha_refs_key)
        alphas_new, v_ratios_nonfi = np.array(alphas_new), np.array(v_ratios_nonfi)


        # calculate the relative error of the alpha_iso
        iso, alpha_aniso = ts_lite(atnames, coords, None, v_ratios_nonfi, alphas_new, beta=beta) 
        rue = np.abs(iso - dft_iso) / dft_iso

        # add to total errors
        tot_rue +=rue

    # return mean relative (absolute) error
    return tot_rue / len(clusters)



if __name__ == "__main__":
    # Tkatchenko Scheffler 2012 for PBE TS-SCS
    # beta = 2.56
    # Ambrosetti 2014 for PBE0 MBD@rsSCS
    beta = 0.85
    lattice = None

    dft_isos_dict = {'BH3': 17.60804, 'BHF2': 16.56286, 'BN': 22.46983, 'BO': 16.82123, 'BS': 34.4018, 'C2H2': 23.44128, 'C2H3': 26.89556, 'C2H4': 28.05114, 'C2H': 24.16172, 'CH2BH': 31.24674, 'CH2F': 16.51907, 'CH2NH': 22.93359, 'CH2PH': 41.66976, 'CH2-t': 15.23122, 'CH3BH2': 29.70554, 'CH3Cl': 30.44706, 'CH3F': 16.89842, 'CH3NH2': 25.85299, 'CH3OH': 21.45095, 'CH3O': 60.7563, 'CH3SH': 37.36721, 'CH3': 16.34304, 'CH4': 16.84694, 'Cl2': 30.79959, 'ClCN': 30.80331, 'ClF': 18.45117, 'CN': 20.92929, 'CO2': 17.12988, 'CO': 13.09801, 'CSO': 34.0532, 'CS': 28.51627, 'F2': 8.6217, 'FCN': 18.48071, 'FCO': 17.65484, 'FNO': 16.99797, 'H2CN': 19.99494, 'H2O': 9.74364, 'H2': 5.48916, 'HBO': 16.73385, 'HBS': 34.35105, 'HCCCl': 37.54696, 'HCCF': 24.02206, 'HCHO': 18.01506, 'HCHS': 34.6727, 'HCl': 17.62394, 'HCN': 17.10007, 'HCONH2': 27.65972, 'HCOOH': 22.5812, 'HCO': 17.14108, 'HCP': 36.21149, 'He': 1.49334, 'HF': 5.72309, 'HNC': 18.31108, 'HNO': 15.35118, 'HNS': 29.89011, 'HO2': 13.44279, 'HOCl': 22.71241, 'HOF': 11.74709, 'HOOH': 15.15929, 'H': 5.20348, 'Li2': 199.80629, 'LiBH4': 31.00643, 'LiCl': 27.04554, 'LiCN': 24.48947, 'Li': 145.05014, 'Mg2': 166.12611, 'Mg': 75.09595, 'N2H2': 19.42433, 'N2H4': 22.921, 'N2': 11.87995, 'Na2': 247.96428, 'NaCl': 32.13845, 'NaCN': 27.95901, 'NaH': 43.03819, 'NaLi': 223.62042, 'Na': 168.29099, 'NCO': 21.05792, 'Ne': 2.81014, 'NH2Cl': 27.10119, 'NH2F': 15.00191, 'NH2OH': 19.08241, 'NH2': 12.31692, 'NH3O': 21.15892, 'NH3': 14.35241, 'NH': 10.1323, 'NOCl': 33.73603, 'NO': 11.69064, 'NP': 28.30619, 'N': 7.79387, 'O2': 10.69618, 'O3': 19.19798, 'OCl2': 38.15706, 'OCl': 20.76964, 'OF2': 14.53111, 'OF': 9.71451, 'OH': 7.63781, 'P2H4': 58.23055, 'P2': 50.02577, 'PH2OH': 33.72788, 'PH2': 29.91048, 'PH3O': 29.88189, 'PH3': 30.93393, 'PH': 28.28455, 'PS': 45.38103, 'P': 26.37305, 'S2H2': 44.64958, 'S2': 40.72003, 'SCl2': 51.61693, 'SCl': 35.81186, 'SF2': 23.98944, 'SF': 22.25586, 'SH2': 24.85501, 'SH': 22.62708, 'SiH3Cl': 42.56196, 'SiH3F': 28.61593, 'SiH3': 34.61421, 'SiH4': 31.91881, 'SiO': 29.62807, 'SO-trip': 23.35157, 'FH-OH': 13.50044, 'H2O-Li': 234.83737, 'AlF': 39.92468, 'Ar': 11.40389, 'BeH2': 20.63737, 'BeH': 34.2047, 'Be': 42.94245, 'BF': 20.38815, 'BH2Cl': 30.06453, 'BH2F': 17.05495, 'BH2': 20.69122, 'SO2': 25.21068, 'LiH': 29.39449}    
    
    clusters_all = []
    # loop through all materials in chargemol.txt and format as cluster dict
    with open('chargemol.txt', 'r') as f:
        clusters_raw = f.read().strip().split('--')[1:]
        for cluster_raw in clusters_raw:
            cluster = {}
            atoms = cluster_raw.strip().split('\n')
            name = atoms[0].split()[0]
            cluster['name'] = name
            cluster['atnums'] = np.array([int(at.split()[2]) for at in atoms])
            cluster['atnames'] = np.array([at.split()[1] for at in atoms])
            cluster['coords'] = ang_to_bohr * np.array([[float(coord) for coord in at.split()[3:6]] for at in atoms])
            cluster['r3s'] = np.array([float(at.split()[-1]) for at in atoms])
            if name != "CH3O":
                clusters_all.append(cluster)
    
    # narrow to only molecules with common elements (must appear at least 5 times in the rest of the set)
    clusters, dft_isos = [], []
    with open("opt_list.txt", "r") as f:
        f.readline()
        for line in f.readlines():
            cluster_name = line.strip()
            for cluster in clusters_all:
                if cluster['name'] == cluster_name:
                    clusters.append(cluster)
                    dft_isos.append(dft_isos_dict[cluster_name])
    print(f"{len(clusters)} molecules included in the optimization")

    # first, set the order of the alpha_refs (the mapping to atnum)
    alpha_refs_key = set()
    for cluster in clusters:
        # get all atnums from all clusters
        alpha_refs_key = alpha_refs_key.union(set(cluster['atnums']))
    alpha_refs_key = list(alpha_refs_key)
    
    # get the initial guess for alpha_refs
    alpha_refs = [get_ref_alpha(atnum) for atnum in alpha_refs_key]
    print("1. Generating initial guess for alphas:")
    print("    Atnum: Alpha")
    for a, atnum in zip(alpha_refs, alpha_refs_key):
        print(f"    {atnum}: {a} Bohr^3")

    alpha_refs_plus_beta = alpha_refs + [beta]

    print("2. Computing initial MRUE:")
    mrue_init = objective_mre_plus_beta(alpha_refs_plus_beta, clusters, dft_isos, alpha_refs_key)
    mrse_init = objective_mrse_plus_beta(alpha_refs_plus_beta, clusters, dft_isos, alpha_refs_key)
    print(f"    Initial Mean Relative Unsigned Error: {mrue_init*100:0.2f}%")
    print(f"    Initial Mean Relative Signed Error: {mrse_init*100:0.2f}%")

    # get optimial alpha_refs by minimizing mean relative unsigned error of isotropic polarizability of clusters
    res = minimize(objective_mre_plus_beta, alpha_refs_plus_beta, method='Nelder-Mead', args=(clusters, dft_isos, alpha_refs_key)) 
    alpha_refs_plus_beta_opt = res.x
    alpha_refs_opt, beta_opt = alpha_refs_plus_beta_opt[:-1], alpha_refs_plus_beta_opt[-1]
    mrue_opt = res.fun
    mrue_opt_check = objective_mre_plus_beta(alpha_refs_plus_beta_opt, clusters, dft_isos, alpha_refs_key)
    mrse_opt = objective_mrse_plus_beta(alpha_refs_plus_beta_opt, clusters, dft_isos, alpha_refs_key)
    print("3. Optimized alphas:")
    print("    Atnum: Alpha")
    for a, atnum in zip(alpha_refs_opt, alpha_refs_key):
        print(f"    {atnum}: {a} Bohr^3")
    print(f"Optimized beta: {beta_opt}")
    print(f"  - MRUE with Optimized alphas: {mrue_opt*100:0.2f}%\n\n")
    print(f"  - MRUE Check with Optimized alphas: {mrue_opt_check*100:0.2f}%\n\n")
    print(f"  - MRSE with Optimized alphas: {mrse_opt*100:0.2f}%\n\n")

    # now report all alphas with optimized reference pols
    with open('ts_opt_analysis_beta_opt.csv', 'w+') as results:
        results.write("Cluster,Iso Opt,XX Opt,XY Opt,XZ Opt,YY Opt,YZ Opt,ZZ Opt\n")
        for cluster in clusters: 
            # non-FI starting guess
            alphas_opt, v_ratios_nonfi = get_all_aim_alpha0(cluster['atnums'], cluster['r3s'], alpha_refs_opt, alpha_refs_key=alpha_refs_key)
            alphas_opt, v_ratios_nonfi = np.array(alphas_opt), np.array(v_ratios_nonfi)
             
            name = cluster['name']
            print(f"Performing TS Analysis on Material {name}\n{'-'*40}\n")
            
            # get TS ad TS-SCS from Manz alphas
            print(f"*** Using the Re-optimized reference polarizabilities\n    {'-'*45}\n")
            iso_pol_ts, iso_pol_scs, alpha_rsscs_mol = ts_analysis(cluster['atnames'], cluster['coords'], lattice, v_ratios_nonfi, alpha0=alphas_opt, beta=beta_opt)
            
            results.write(f"{name},{iso_pol_scs},{alpha_rsscs_mol[0,0]},{alpha_rsscs_mol[0,1]},{alpha_rsscs_mol[0,2]},{alpha_rsscs_mol[1,1]},{alpha_rsscs_mol[1,2]},{alpha_rsscs_mol[2,2]}\n")


