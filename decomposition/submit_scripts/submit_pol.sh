#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=0
#SBATCH --ntasks-per-node=28
#SBATCH --time=0-11:55 # time (DD-HH:MM)
#SBATCH --job-name=MgO_cluster # Name of job in queue
#SBATCH --account=rrg-ipaci # Paci group resource allocation

module load StdEnv/2020  gcc/10.3.0  openmpi/4.1.1
module load orca/5.0.4L

job=${SLURM_JOB_NAME}
job=$(echo ${job%%.*})

#Start ORCA job. ORCA is started using full pathname (necessary for parallel execution). Output file is written directly to submit directory on frontnode.
echo "Running Zero Field"
${EBROOTORCA}/orca input.inp >> output.out

echo "Running Xfield"
cd xfield
${EBROOTORCA}/orca input.inp >> output.out
echo "Running Yfield"
cd ../yfield
${EBROOTORCA}/orca input.inp >> output.out
echo "Running Zfield"
cd ../zfield
${EBROOTORCA}/orca input.inp >> output.out
echo "Running Xfield"

echo "Running Minus Xfield"
cd ../x-field
${EBROOTORCA}/orca input.inp >> output.out
echo "Running Minus Yfield"
cd ../y-field
${EBROOTORCA}/orca input.inp >> output.out
echo "Running Minus Zfield"
cd ../z-field
${EBROOTORCA}/orca input.inp >> output.out

