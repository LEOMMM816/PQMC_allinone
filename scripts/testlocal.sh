#!/bin/bash
# this script is used to run the code locally for testing
# it will create a work directory, compile the code, run it with mpirun,
# and then post-process the results
# this file is run from the scripts/ directory

# ------------------------

# delete the workplace directory if it exists
rm -rf ../workplace

# create the work directory and copy everything needed
submit_dir=$(pwd)/..
work_dir=$submit_dir/workplace
mkdir $work_dir
mkdir -p $work_dir/data/data $work_dir/data/out_files
mkdir -p $work_dir/build $work_dir/bin

cp -r $submit_dir/src $work_dir/
cp -r $submit_dir/input $work_dir/
cp -r $submit_dir/scripts $work_dir/
cp -r 
# ------------------------
: "${NBLOCK:?NBLOCK not set}"   # 如果没传 NBLOCK 就报错退出
echo "In test.sh: NBLOCK = $NBLOCK"
nperblock=1
ntasks=$(( NBLOCK * nperblock ))

# run the fortran code
cd ./code
rm -rf ./*.mod
rm -rf ./*.out
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
mpifort mod_nrtype.f90 mod_nrutil.f90  mod_matrixlib.f90  mod_ranstate.f90 input.f90 mod_lattice.f90 mod_phonon_field.f90 mod_evolution.f90 mod_update.f90 mod_meas.f90 Main_PQMC.f90 \
 -o main.out -cpp -DMPI -DCPPBLOCK=$NBLOCK -DCPPPERBLOCK=$nperblock -lopenblas -g -march=native -fno-omit-frame-pointer -O3 -fbacktrace -ffpe-trap=invalid,zero,overflow -Wno-unused-parameter -Wno-unused-variable -finit-real=snan
cp main.out ../
cd ..
 mpirun -np "$ntasks" ./script/run_perf.sh ./main.out -> main.log
#mpirun -np "$ntasks" ./main.out > main.log
rm -rf ./namelists

# post process
cd ./code
gfortran mod_nrtype.f90 outputnew.f90 -cpp -DMPI -fcheck=all -g -o pp.out
cp pp.out ../
cd ..
./pp.out > pp.log

