#!/bin/bash
#SBATCH -p intel-sc3-32c
#SBATCH -q normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH -J ctL40Oi
#SBATCH -x intel124
#SBATCH --mem-per-cpu=5G
hostname
ntasks=32
jobname=ctL40Oi
output_dir=$(pwd)/results
dir_local=/data/$$               #tmp directory private to specific excute compute node  
mkdir -p $dir_local              #$$ refer to  unique pid of spccific job

module load gcc/9.4.0
module load openmpi/4.1.0
export OPENBLAS_NUM_THREADS=1 
module load OpenBLAS/0.3.13 
mpifort input.f90 mod_matrixlib.f90 mod_nrtype.f90 mod_nrutil.f90 mod_ranstate.f90 mod_lattice.f90 mod_phonon_field.f90 mod_evolution.f90 mod_update.f90 mod_gpt_meas.f90 Main_PQMC.f90 \
 -cpp -DMPI -lopenblas -g -O3 
cp ./Creutz_holstein.nml $dir_local
cp ./a.out $dir_local
cp ./cp_middle.sh $dir_local
cp ./input.f90 $dir_local
cp ./outputfinal.f90 $dir_local
cd $dir_local #now we are utilizing the SSD disk locate at compute node
export OPENBLAS_NUM_THREADS=1

# ------------------------
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
mkdir -p ./data/data ./data/middata ./data/ML_data ./data/raw_data
mpirun -np 32 ./a.out > out.log 

# go to the output directory and output final results

#gfortran input.f90 outputfinal.f90 -cpp -DMPI 
#./a.out $jobname

# copy file and remove tmp directory on node
rm -r namelists
cd /data
cp -r $dir_local $output_dir
rm -rf $dir_local
