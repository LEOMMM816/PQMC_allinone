#!/bin/bash

ntasks=${SLURM_TASKS_PER_NODE}
nperblock=$(( ntasks / NBLOCK ))
jobname=${SLURM_JOB_NAME}
output_dir=$(pwd)/results/${jobname}
dir_local=/data/$$               #tmp directory private to specific excute compute node  
mkdir -p $dir_local              #$$ refer to  unique pid of spccific job
module load gcc/9.4.0
module load openmpi/4.1.0
module load OpenBLAS/0.3.13 
export OPENBLAS_NUM_THREADS=1
# delete the data directory if it exists
rm -rf ./data
rm -rf ./*.out
rm -rf ./*.log
rm -rf ./results
# create the data directory

# ------------------------
# compile the fortran code
cd ./code
mpifort input.f90 mod_nrtype.f90 mod_nrutil.f90 mod_matrixlib.f90  mod_ranstate.f90 mod_lattice.f90 mod_phonon_field.f90 mod_evolution.f90 mod_update.f90 mod_meas.f90 Main_PQMC.f90 \
 -cpp -DMPI -DCPPBLOCK=$NBLOCK -DCPPPERBLOCK=$nperblock  -lopenblas -g -O3 -o main.out
cp ./main.out ${dir_local}/
cd ..
# copy necessary files to the tmp directory on node
cp -r ./script/cp_middle.sh ${dir_local}/
cp -r ./code ${dir_local}/
mkdir -p ${dir_local}/model
cp -r ./model/temp ${dir_local}/model/
cd $dir_local #now we are utilizing the SSD disk locate at compute node
mkdir -p ./data/data ./data/out_files
export OPENBLAS_NUM_THREADS=1

# run the fortran code
mpirun -np "$ntasks" ./main.out > main.log 

# post process
cd ./code
gfortran outputnew.f90 -cpp -DMPI -fcheck=all -g -o pp.out
cp pp.out ../
cd ..
./pp.out > pp.log
# copy the whole directory back to storage directory
cd /data
mkdir -p $output_dir
cp -a $dir_local $output_dir
rm -rf $dir_local
