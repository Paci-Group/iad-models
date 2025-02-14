import numpy as np

BOHR = 1.88972598858
VALENCE = {"O": 8, "Mg": 12, "Ti": 22, "Si": 14, "H": 1, "Ag": 19}

def get_valence_from_xyz(xyz):
    valence = []
    with open(xyz, 'r') as f:
        nat = int(f.readline().strip())
        f.readline()
        for _ in range(nat):
            valence.append(VALENCE[f.readline().strip().split()[0]])
    return np.array(valence)

def parse_analytical_pol(outfile):
    raw_tensor = np.zeros((3,3))
    diagonalized_tensor = np.zeros(3)
    isotropic_pol = 0
    with open(outfile, 'r') as f:
        line = f.readline()
        while line:
            if line.startswith("The raw cartesian tensor"):
                for i in range(3):
                    spline = f.readline().strip().split()
                    raw_tensor[i, :] = np.array([float(e) for e in spline])
            elif line.startswith("diagonalized tensor"):
                diagonalized_tensor[:] = np.array([float(e) for e in f.readline().strip().split()])
            elif line.startswith("Isotropic polarizability"):
                isotropic_pol = float(line.strip().split()[3])
            line = f.readline()
    return raw_tensor, diagonalized_tensor, isotropic_pol

def parse_dipole(outfile):
    with open(outfile, 'r') as f:
        for line in f.readlines():
            if line.startswith("Total Dipole Moment"):
                d = np.array([float(e) for e in line.strip().split()[4:]])
    return d

def parse_numerical_pol(outfile_f1, outfile_f2, ef=0.001):
    d1 = parse_dipole(outfile_f1)
    d2 = parse_dipole(outfile_f2)
    return (d2 - d1) / ef

def parse_moments(outfile, core_charges):
    moments = {"positions": None, "charges": None, "dipoles": None}
    with open(outfile, 'r') as f:
        line = f.readline()
        while line:
            if line.startswith("  Number of atoms:"):
                nat = int(line.strip().split()[-1])
            elif line.startswith("# at         x                y                z"):
                pos = []
                for _ in range(nat):
                    line = f.readline()
                    pos.append([float(e) for e in line.strip().split()[1:4]])
                moments["positions"] = np.array(pos) * BOHR
            elif line.startswith("# Id   cp   ncp   Name  Z           1               x               z               y"):
                charges = np.zeros(nat)
                dipoles = np.zeros((nat,3))
                for i in range(nat):
                    line = f.readline()
                    spline = line.strip().split()
                    charges[i] = float(spline[5])
                    dipoles[i, :] = np.array([float(spline[6]), float(spline[8]), float(spline[7])])
                moments["charges"] = core_charges - charges
                moments["dipoles"] = -dipoles
            line = f.readline()
    return moments

def pretty_print_cluster_pols(raw_ana, diag_ana, iso_ana, num_ex, num_ey, num_ez, num_iso, f=None):
    print(f"* Analyzing Cluster {cluster}\n{'-'*50}", file=f)
    print("  - Extracting cluster polarizabilities from analytical and finite difference calculations.\n", file=f)
    print(f"    Analytical Polarizability Tensor\n    {'-'*50}", file=f)
    for row in raw_ana:
        print(f"    {row[0]:.5f}    {row[1]:.5f}    {row[2]:.5f}", file=f)
    print(f"\n    Diagonalized Analytical Polarizability\n    {'-'*50}", file=f)
    print(f"    {diag_ana[0]:.5f}    {diag_ana[1]:.5f}    {diag_ana[2]:.5f}\n", file=f)
    print(f"    Isotropic Analytical Polarizability: {iso_ana:.5f}\n", file=f) 
    
    print(f"    Numerical Polarizability Tensor\n    {'-'*50}", file=f)
    raw_num = np.concatenate([num_ex.reshape(-1,1), num_ey.reshape(-1,1), num_ez.reshape(-1,1)], axis=1)
    for row in raw_num:
        print(f"    {row[0]:.5f}    {row[1]:.5f}    {row[2]:.5f}", file=f)

    diag_num, _ = np.linalg.eig(raw_num)
    diag_num = np.sort(diag_num)
    print(f"\n    Diagonalized Numerical Polarizability\n    {'-'*50}", file=f)
    print(f"    {diag_num[0]:.5f}    {diag_num[1]:.5f}    {diag_num[2]:.5f}\n", file=f)
    print(f"    Isotropic Numerical Polarizability: {num_iso:.5f}\n", file=f) 


def pretty_print_atomic_pols(positions, d_mup_dex, d_mup_dey, d_mup_dez, d_muq_dex, d_muq_dey, d_muq_dez, f=None):
    print("  - Extracting site-specific polarizabilities from Bader analysis.\n", file=f)
    print("    +++ reporting only the diagonal elements of the polarizability tensors for each atom", file=f)
    print("    * Atomic Positions (in Bohr)", file=f)
    for pos in positions:
        print(f"      {pos[0]:.5f}    {pos[1]:.5f}    {pos[2]:.5f}", file=f)
    print("    * Local Dipole Polarizabilities (Bohr ^3)", file=f)
    for apx, apy, apz in zip(d_mup_dex, d_mup_dey, d_mup_dez):
        print(f"      {apx[0]:.5f}    {apy[1]:.5f}    {apz[2]:.5f}", file=f)
    print("    * Charge Transfer Polarizabilities (Bohr ^3)", file=f)
    for aqx, aqy, aqz in zip(d_muq_dex, d_muq_dey, d_muq_dez):
        print(f"      {aqx[0]:.5f}    {aqy[1]:.5f}    {aqz[2]:.5f}", file=f)
    print("    * Total Polarizabilities (Bohr ^3)", file=f)
    for apx, apy, apz, aqx, aqy, aqz in zip(d_mup_dex, d_mup_dey, d_mup_dez, d_muq_dex, d_muq_dey, d_muq_dez):
        print(f"      {apx[0]+aqx[0]:.5f}    {apy[1]+aqy[1]:.5f}    {apz[2]+aqz[2]:.5f}", file=f)


def pretty_print_decomposed_cluster_pols(alpha_p_x, alpha_p_y, alpha_p_z, alpha_q_x, alpha_q_y, alpha_q_z, f=None):
    print("\n    +++ Summing over atomic polarizabilities to obtain decomposed cluster polarizabilities.", file=f)
    iso_alpha_p = (alpha_p_x[0] + alpha_p_y[1] + alpha_p_z[2]) / 3 
    iso_alpha_q = (alpha_q_x[0] + alpha_q_y[1] + alpha_q_z[2]) / 3
    
    print(f"    Summed Atomic Dipole Polarizability Tensor\n    {'-'*50}", file=f)
    raw_p = np.concatenate([alpha_p_x.reshape(-1,1), alpha_p_y.reshape(-1,1), alpha_p_z.reshape(-1,1)], axis=1)
    for row in raw_p:
        print(f"    {row[0]:.5f}    {row[1]:.5f}    {row[2]:.5f}", file=f)
    diag_p, _ = np.linalg.eig(raw_p)
    diag_p = np.sort(diag_p)
    print(f"\n    Diagonalized Summed Atomic Dipole Polarizability\n    {'-'*50}", file=f)
    print(f"    {diag_p[0]:.5f}    {diag_p[1]:.5f}    {diag_p[2]:.5f}\n", file=f)

    print(f"    Summed Atomic Charge Transfer Polarizability Tensor\n    {'-'*50}", file=f)
    raw_q = np.concatenate([alpha_q_x.reshape(-1,1), alpha_q_y.reshape(-1,1), alpha_q_z.reshape(-1,1)], axis=1)
    for row in raw_q:
        print(f"    {row[0]:.5f}    {row[1]:.5f}    {row[2]:.5f}", file=f)
    diag_q, _ = np.linalg.eig(raw_q)
    diag_q = np.sort(diag_q)
    print(f"\n    Diagonalized Summed Atomic Charge Transfer Polarizability\n    {'-'*50}", file=f)
    print(f"    {diag_q[0]:.5f}    {diag_q[1]:.5f}    {diag_q[2]:.5f}\n", file=f)
    
    print(f"    Total Summed Atomic Tensor\n    {'-'*50}", file=f)
    for row in raw_p+raw_q:
        print(f"    {row[0]:.5f}    {row[1]:.5f}    {row[2]:.5f}", file=f)
    diag_pq, _ = np.linalg.eig(raw_p + raw_q)
    diag_pq = np.sort(diag_pq)
    print(f"\n    Diagonalized Total Summed Atomic Polarizability\n    {'-'*50}", file=f)
    print(f"    {diag_pq[0]:.5f}    {diag_pq[1]:.5f}    {diag_pq[2]:.5f}\n", file=f)

    print(f"    Isotropic Summed Atomic Dipole Polarizability: {iso_alpha_p:.5f}", file=f) 
    print(f"    Isotropic Summed Atomic Charge Transfer Polarizability: {iso_alpha_q:.5f}", file=f) 
    print(f"    Isotropic Total Summed Atomic Polarizability: {iso_alpha_p+iso_alpha_q:.5f}", file=f) 



if __name__ == "__main__":
    ef = 0.001 * 2
    clusters = [f"tio2_{n}" for n in range(2,11)]
    with open("zero_centered_site_specific_analysis.txt", "w+") as f:
        print("Cluster,XX Analytical,XY Analytical,XZ Analytical,YY Analytical,YZ Analytical,ZZ Analytical,Iso Analytical,D1 Analytical,D2 Analytical,D3 Analytical,XX Numerical,XY Numerical,XZ Numerical,YY Numerical,YZ Numerical,ZZ Numerical,Iso Numerical,D1 Numerical,D2 Numerical,D3 Numerical,XX P,XY P,XZ P,YY P,YZ P,ZZ P,XX Q,XY Q,XZ Q,YY Q,YZ Q,ZZ Q,XX P+Q,XY P+Q,XZ P+Q,YY P+Q,YZ P+Q,ZZ P+Q,Iso P,Iso Q,Iso P+Q")
        for cluster in clusters:
            valence = get_valence_from_xyz(f"{cluster}/input.xyz")
            raw_ana, diag_ana, iso_ana = parse_analytical_pol(f"{cluster}/output.out")

            num_ex = parse_numerical_pol(f"{cluster}/x-field/output.out", f"{cluster}/xfield/output.out",ef=ef) 
            num_ey = parse_numerical_pol(f"{cluster}/y-field/output.out", f"{cluster}/yfield/output.out",ef=ef) 
            num_ez = parse_numerical_pol(f"{cluster}/z-field/output.out", f"{cluster}/zfield/output.out",ef=ef) 
            num_iso = (num_ex[0] + num_ey[1] + num_ez[2]) / 3
            raw_num = np.concatenate([num_ex.reshape(-1,1), num_ey.reshape(-1,1), num_ez.reshape(-1,1)], axis=1)
            diag_num, _ = np.linalg.eig(raw_num) 
            diag_num = np.sort(diag_num)

            #print information about overal cluster polarizability
            pretty_print_cluster_pols(raw_ana, diag_ana, iso_ana, num_ex, num_ey, num_ez, num_iso, f=f) 

            moments_nex = parse_moments(f"{cluster}/x-field/output.cro", valence)
            moments_ney = parse_moments(f"{cluster}/y-field/output.cro", valence)
            moments_nez = parse_moments(f"{cluster}/z-field/output.cro", valence)
            moments_ex = parse_moments(f"{cluster}/xfield/output.cro", valence)
            moments_ey = parse_moments(f"{cluster}/yfield/output.cro", valence)
            moments_ez = parse_moments(f"{cluster}/zfield/output.cro", valence)
            com = moments_ex["positions"].mean(axis=0)

            # local dipole moment polarizability
            d_mup_dex = (moments_ex["dipoles"] - moments_nex["dipoles"]) / ef
            d_mup_dey = (moments_ey["dipoles"] - moments_ney["dipoles"]) / ef
            d_mup_dez = (moments_ez["dipoles"] - moments_nez["dipoles"]) / ef
            
            # charge dipole moment
            muq_nex = moments_nex["charges"].reshape(-1,1) * (moments_nex["positions"] - com)
            muq_ney = moments_ney["charges"].reshape(-1,1) * (moments_ney["positions"] - com)
            muq_nez = moments_nez["charges"].reshape(-1,1) * (moments_nez["positions"] - com)
            muq_ex = moments_ex["charges"].reshape(-1,1) * (moments_ex["positions"] - com)
            muq_ey = moments_ey["charges"].reshape(-1,1) * (moments_ey["positions"] - com)
            muq_ez = moments_ez["charges"].reshape(-1,1) * (moments_ez["positions"] - com)
            
            # charge transfer polarizability
            d_muq_dex = (muq_ex - muq_nex) / ef
            d_muq_dey = (muq_ey - muq_ney) / ef
            d_muq_dez = (muq_ez - muq_nez) / ef

            pretty_print_atomic_pols(moments_ex["positions"], d_mup_dex, d_mup_dey, d_mup_dez, d_muq_dex, d_muq_dey, d_muq_dez, f=f)

            # validation: sum of atomic polarizabilities should equal total cluster polarizability
            alpha_p_x = d_mup_dex.sum(axis=0)
            alpha_p_y = d_mup_dey.sum(axis=0)
            alpha_p_z = d_mup_dez.sum(axis=0)
            
            alpha_q_x = d_muq_dex.sum(axis=0)
            alpha_q_y = d_muq_dey.sum(axis=0)
            alpha_q_z = d_muq_dez.sum(axis=0)
            
            iso_alpha_p = (alpha_p_x[0] + alpha_p_y[1] + alpha_p_z[2]) / 3 
            iso_alpha_q = (alpha_q_x[0] + alpha_q_y[1] + alpha_q_z[2]) / 3
            
            pretty_print_decomposed_cluster_pols(alpha_p_x, alpha_p_y, alpha_p_z, alpha_q_x, alpha_q_y, alpha_q_z, f=f)
            print("\n\n", file=f)

            print(f"{cluster},{raw_ana[0,0]:.5f},{raw_ana[0,1]:.5f},{raw_ana[0,2]:.5f},{raw_ana[1,1]:.5f},{raw_ana[1,2]:.5f},{raw_ana[2,2]:.5f},{iso_ana:.5f},{diag_ana[0]:.5f},{diag_ana[1]:.5f},{diag_ana[2]:.5f},{num_ex[0]:.5f},{(num_ex[1]+num_ey[0])/2:.5f},{(num_ex[2]+num_ez[0])/2:.5f},{num_ey[1]:.5f},{(num_ey[2]+num_ez[1])/2:.5f},{num_ez[2]:.5f},{num_iso:.5f},{diag_num[0]:.5f},{diag_num[1]:.5f},{diag_num[2]:.5f},{alpha_p_x[0]:.5f},{(alpha_p_x[1]+alpha_p_y[0])/2:.5f},{(alpha_p_x[2]+alpha_p_z[0])/2:.5f},{alpha_p_y[1]:.5f},{(alpha_p_y[2]+alpha_p_z[1])/2:.5f},{alpha_p_z[2]:.5f},{alpha_q_x[0]:.5f},{(alpha_q_x[1]+alpha_q_y[0])/2:.5f},{(alpha_q_x[2]+alpha_q_z[0])/2:.5f},{alpha_q_y[1]:.5f},{(alpha_q_y[2]+alpha_q_z[1])/2:.5f},{alpha_q_z[2]:.5f},{alpha_p_x[0]+alpha_q_x[0]:.5f},{(alpha_p_x[1]+alpha_p_y[0]+alpha_q_x[1]+alpha_q_y[0])/2:.5f},{(alpha_p_x[2]+alpha_p_z[0]+alpha_q_x[2]+alpha_q_z[0])/2:.5f},{alpha_p_y[1]+alpha_q_y[1]:.5f},{(alpha_p_y[2]+alpha_p_z[1]+alpha_q_y[2]+alpha_q_z[1])/2:.5f},{alpha_p_z[2]+alpha_q_z[2]:.5f},{iso_alpha_p:.5f},{iso_alpha_q:.5f},{iso_alpha_p+iso_alpha_q:.5f}")

