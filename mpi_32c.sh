#!/bin/bash
#SBATCH -p intel-sc3-32c
#SBATCH -q normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH -J ctL40O2cp
#SBATCH -x intel124
#SBATCH --mem-per-cpu=5G
hostname
jobname=ctL40O2cp
output_dir=$(pwd)/results
dir_local=/data/$$               #tmp directory private to specific excute compute node  
mkdir -p $dir_local              #$$ refer to  unique pid of spccific job

module load gcc/9.3.0
module load openmpi/4.1.0
export OPENBLAS_NUM_THREADS=1 
module load OpenBLAS/0.3.13 
mpifort input.f90 mod_matrixlib.f90 mod_nrtype.f90 mod_nrutil.f90 mod_ranstate.f90 mod_zxy_pqmc_ssh_server.f90 mod_pqmc_zxy_measure.f90 hubburd_pqmc_ssh_server.f90 -O3 -DMPI -cpp -lopenblas -fcheck=all -g
cp ./a.out $dir_local
cp ./cp_middle.sh $dir_local
cp ./input.f90 $dir_local
cp ./outputfinal.f90 $dir_local
cd $dir_local #now we are utilizing the SSD disk locate at compute node
export OPENBLAS_NUM_THREADS=1
mkdir -p ./data/data ./data/middata ./data/ML_data ./data/raw_data
mpirun -np 32 ./a.out > out.log 

# go to the output directory and output final results

gfortran input.f90 outputfinal.f90 -cpp -DMPI 
./a.out $jobname

# copy file and remove tmp directory on node
cd /data
cp -r $dir_local $output_dir
rm -rf $dir_local
