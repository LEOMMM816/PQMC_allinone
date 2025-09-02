#!/bin/bash
# delete the data directory if it exists
rm -rf ./data
# create the data directory
mkdir -p ./data/data ./data/middata ./data/ML_data ./data/raw_data
# run the fortran code
mpifort input.f90 mod_matrixlib.f90 mod_nrtype.f90 mod_nrutil.f90 mod_ranstate.f90 mod_lattice.f90 mod_phonon_field.f90 mod_evolution.f90 mod_update.f90 mod_gpt_meas.f90 Main_PQMC.f90 \
 -cpp -DMPI -lopenblas -g -O0 -fbacktrace -ffpe-trap=invalid,zero,overflow -Wall -Wextra -Wno-unused-parameter -Wno-unused-variable -finit-real=snan
mpirun -np 1 -quiet ./a.out

#gfortran input.f90 outputfinal.f90 -cpp -DMPI -fcheck=all -g
#./a.out test > out.log

