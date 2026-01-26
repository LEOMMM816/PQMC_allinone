#!/bin/bash

echo "----- post processing -----"
cd ../code
gfortran mod_nrtype.f90 outputnew.f90 -cpp -DMPI -fcheck=all -g -o pp.out || true
cp pp.out ../ 
cd ..
./pp.out &> pp.log 
cd ./script