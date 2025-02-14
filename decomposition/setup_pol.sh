#!/bin/bash

for g in geometries/*.xyz
do
  echo $g
  name=${g##*/}
  name="${name%.xyz}"
  echo "Setting up Job for $name"
  
  # zero field
  mkdir $name
  cp $g $name/input.xyz
  cp template_zerofield.inp $name/input.inp
  
  # plus field
  mkdir $name/zfield $name/yfield $name/xfield
  cp $name/input.xyz $name/xfield
  cp $name/input.xyz $name/yfield
  cp $name/input.xyz $name/zfield
  cp template_xfield.inp $name/xfield/input.inp
  cp template_yfield.inp $name/yfield/input.inp
  cp template_zfield.inp $name/zfield/input.inp
  
  # minus field
  mkdir $name/z-field $name/y-field $name/x-field
  cp $name/input.xyz $name/x-field
  cp $name/input.xyz $name/y-field
  cp $name/input.xyz $name/z-field
  cp $tmps/template_x-field.inp $name/x-field/input.inp
  cp $tmps/template_y-field.inp $name/y-field/input.inp
  cp $tmps/template_z-field.inp $name/z-field/input.inp
  
  # submit
  cp submit_pol.sh $name/submit.sh
  cd $name
  sbatch submit.sh
  cd ..
done
