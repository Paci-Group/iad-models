#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=0
#SBATCH --ntasks-per-node=28
#SBATCH --time=0-11:55 # time (DD-HH:MM)
#SBATCH --job-name=tio2_critic2 # Name of job in queue
#SBATCH --account=rrg-ipaci # Paci group resource allocation

job=${SLURM_JOB_NAME}
job=$(echo ${job%%.*})
tmps="templates"

for d in tio2_*
do 
  cp $tmps/critic_input_1.cri $d/input.cri
  cp $tmps/critic_input_1.cri $d/xfield/input.cri
  cp $tmps/critic_input_1.cri $d/yfield/input.cri
  cp $tmps/critic_input_1.cri $d/zfield/input.cri
  #cp $tmps/critic_input_2.cri $d/xfield/input.neutbasin.cri
  #cp $tmps/critic_input_2.cri $d/yfield/input.neutbasin.cri
  #cp $tmps/critic_input_2.cri $d/zfield/input.neutbasin.cri
  cp $tmps/critic_input_1.cri $d/x-field/input.cri
  cp $tmps/critic_input_1.cri $d/y-field/input.cri
  cp $tmps/critic_input_1.cri $d/z-field/input.cri
  #cp $tmps/critic_input_2.cri $d/x-field/input.neutbasin.cri
  #cp $tmps/critic_input_2.cri $d/y-field/input.neutbasin.cri
  #cp $tmps/critic_input_2.cri $d/z-field/input.neutbasin.cri
  cd $d
  echo "Running Critic2 for $d"
  critic2 input.cri > output.cro
  echo "  - Running Critic2 for $d/xfield"
  cd xfield
  critic2 input.cri > output.cro
  #critic2 input.neutbasin.cri > output.neutbasin.cro
  echo "  - Running Critic2 for $d/yfield"
  cd ../yfield
  critic2 input.cri > output.cro
  #critic2 input.neutbasin.cri > output.neutbasin.cro
  echo "  - Running Critic2 for $d/zfield"
  cd ../zfield
  critic2 input.cri > output.cro
  #critic2 input.neutbasin.cri > output.neutbasin.cro
  echo "  - Running Critic2 for $d/x-field"
  cd ../x-field
  critic2 input.cri > output.cro
  #critic2 input.neutbasin.cri > output.neutbasin.cro
  echo "  - Running Critic2 for $d/y-field"
  cd ../y-field
  critic2 input.cri > output.cro
  #critic2 input.neutbasin.cri > output.neutbasin.cro
  echo "  - Running Critic2 for $d/z-field"
  cd ../z-field
  critic2 input.cri > output.cro
  #critic2 input.neutbasin.cri > output.neutbasin.cro
  cd ../../
done
