#!/bin/bash
# delete the data directory if it exists
rm -rf ./data
# create the data directory
mkdir -p ./data/data ./data/middata ./data/ML_data ./data/raw_data
# run the fortran code
mpifort input.f90 mod_matrixlib.f90 mod_nrtype.f90 mod_nrutil.f90 mod_ranstate.f90 mod_zxy_pqmc_ssh_server.f90 mod_pqmc_zxy_measure.f90 hubburd_pqmc_ssh_server.f90 -O3 -DMPI -cpp -lopenblas -fcheck=all -g
mpirun -np 4 ./a.out > out.log

gfortran input.f90 outputfinal.f90 -cpp -DMPI -fcheck=all -g
./a.out test > out.log

