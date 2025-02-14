#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=0
#SBATCH --ntasks-per-node=28
#SBATCH --time=0-11:55 # time (DD-HH:MM)
#SBATCH --job-name=mgo_cubegen # Name of job in queue
#SBATCH --account=rrg-ipaci # Paci group resource allocation

module load StdEnv/2020  gcc/10.3.0  openmpi/4.1.1
module load orca/5.0.4L

job=${SLURM_JOB_NAME}
job=$(echo ${job%%.*})
tmps="templates"

for d in tio2_*
do  
  cp $tmps/get_cube.inp $d/cubegen.inp
  cp $tmps/get_cube.inp $d/xfield/cubegen.inp
  cp $tmps/get_cube.inp $d/yfield/cubegen.inp
  cp $tmps/get_cube.inp $d/zfield/cubegen.inp
  cd $d
  echo "Getting cube for $d"
  ${EBROOTORCA}/orca cubegen.inp >> cubegen.out
  echo "  - Running Critic2 for $d/x-field"
  cd xfield
  ${EBROOTORCA}/orca cubegen.inp >> cubegen.out
  echo "  - Running Critic2 for $d/y-field"
  cd ../yfield
  ${EBROOTORCA}/orca cubegen.inp >> cubegen.out
  echo "  - Running Critic2 for $d/z-field"
  cd ../zfield
  ${EBROOTORCA}/orca cubegen.inp >> cubegen.out
  cd ../../
done
