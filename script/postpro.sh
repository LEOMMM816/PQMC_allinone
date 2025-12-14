#!/bin/bash

echo "----- post processing -----"
cd ./code
gfortran outputnew.f90 -cpp -DMPI -fcheck=all -g -o pp.out || true
cp pp.out ../ || true
cd ..
./pp.out &> pp.log || true