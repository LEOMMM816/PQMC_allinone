#!/bin/bash
# delete the data directory if it exists
rm -rf ./data
# create the data directory
mkdir -p ./data/data ./data/out_files
# ------------------------
ntasks=4
TEMPLATE=Creutz_holstein.nml        # your source file
PREFIX=input            # target prefix
W=4                          # zero-pad width -> 0000..0127
OUTDIR=namelists             # new folder to hold the copies

# make the folder (safe if it already exists)
mkdir -p "$OUTDIR"

# loop and copy
for i in $(seq 0 $((ntasks-1))); do
  printf -v tag "%0${W}d" "$i"      # e.g. 0000, 0001, ...
  cp "$TEMPLATE" "$OUTDIR/${PREFIX}${tag}.nml"
done
#-------------------------
# run the fortran code
mpifort input.f90 mod_matrixlib.f90 mod_nrtype.f90 mod_nrutil.f90 mod_ranstate.f90 mod_lattice.f90 mod_phonon_field.f90 mod_evolution.f90 mod_update.f90 mod_gpt_meas.f90 Main_PQMC.f90 \
 -cpp -DMPI -lopenblas -g -O0 -fbacktrace -ffpe-trap=invalid,zero,overflow -Wall -Wextra -Wno-unused-parameter -Wno-unused-variable -finit-real=snan
mpirun -np 4 ./a.out -> out.log
rm -rf ./namelists
# post process
gfortran outputnew.f90 -cpp -DMPI -fcheck=all -g
./a.out > out.log

