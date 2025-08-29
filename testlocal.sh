#!/bin/bash
# delete the data directory if it exists
rm -rf ./data
# create the data directory
mkdir -p ./data/data ./data/middata ./data/ML_data ./data/raw_data
# run the fortran code
mpifort input.f90 mod_matrixlib.f90 mod_nrtype.f90 mod_nrutil.f90 mod_ranstate.f90 mod_lattice.f90 mod_phonon_field.f90 mod_evolution.f90 mod_update.f90 mod_measure.f90 Main_PQMC.f90 -cpp -DMPI -lopenblas -fcheck=all -g
mpirun -np 1 ./a.out

#gfortran input.f90 outputfinal.f90 -cpp -DMPI -fcheck=all -g
#./a.out test > out.log

