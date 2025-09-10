#!/bin/bash
# delete the data directory if it exists
rm -rf ./data
rm -rf ./*.out
rm -rf ./*.log
# create the data directory
mkdir -p ./data/data ./data/out_files
# ------------------------
nblock=4
nperblock=2
ntasks=$(( nblock * nperblock ))
TEMPLATE=model/Creutz_holstein.nml        # your source file
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
cd ./code
mpifort input.f90 mod_matrixlib.f90 mod_nrtype.f90 mod_nrutil.f90 mod_ranstate.f90 mod_lattice.f90 mod_phonon_field.f90 mod_evolution.f90 mod_update.f90 mod_meas.f90 Main_PQMC.f90 \
 -o main.out -cpp -DMPI -DCPPBLOCK=$nblock -DCPPPERBLOCK=$nperblock -lopenblas -g -O0 -fbacktrace -ffpe-trap=invalid,zero,overflow -Wall -Wextra -Wno-unused-parameter -Wno-unused-variable -finit-real=snan
cp main.out ../
cd ..
mpirun -np "$ntasks" ./main.out -> main.log
rm -rf ./namelists

# post process
cd ./code
gfortran outputnew.f90 -cpp -DMPI -fcheck=all -g -o pp.out
cp pp.out ../
cd ..
./pp.out > pp.log

