#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=0
#SBATCH --ntasks-per-node=28
#SBATCH --time=0-11:55 # time (DD-HH:MM)
#SBATCH --job-name=tio2_chargemol # Name of job in queue
#SBATCH --account=rrg-ipaci # Paci group resource allocation

module load StdEnv/2020  gcc/10.3.0  openmpi/4.1.1
module load orca/5.0.4L

job=${SLURM_JOB_NAME}
job=$(echo ${job%%.*})
tmps="templates"

for d in tio2_*/
do 
  cp $tmps/molden2wfx.in $d/molden2wfx.in
  cd $d
  echo "Working on $d"
  echo "  - Running orca_2mkl"
  ${EBROOTORCA}/orca_2mkl input -molden
  echo "  - Converting molden file to wfx with MultiWFN"
  multiwfn input.molden.input -silent < molden2wfx.in > output.log
  rm output.log
  cp ../$tmps/job_control.txt .
  echo "  - Running Chargemol"
  chargemol_parallel
  rm job_control.txt chargemol.wfx
  rm *molden*
  rm DDEC6_even_tempered_bond_orders.xyz DDEC_atomic_Rfourth_moments.xyz  DDEC_atomic_Rsquared_moments.xyz
  cd ../
done
