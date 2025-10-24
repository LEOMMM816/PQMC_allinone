#!/bin/bash
#SBATCH -p intel-sc3-32c
#SBATCH -q normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH -J epzL10Oi
#SBATCH -x intel124
#SBATCH --mem-per-cpu=5G
hostname
nblock=8
nperblock=8
ntasks=$(( nblock * nperblock ))
jobname=epzL10Oi
output_dir=$(pwd)/results
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
# create the data directory
mkdir -p ./data/data ./data/out_files
# ------------------------
# compile the fortran code
cd ./code
mpifort input.f90 mod_nrtype.f90 mod_nrutil.f90 mod_matrixlib.f90  mod_ranstate.f90 mod_lattice.f90 mod_phonon_field.f90 mod_evolution.f90 mod_update.f90 mod_meas.f90 Main_PQMC.f90 \
 -cpp -DMPI -DCPPBLOCK=$nblock -DCPPPERBLOCK=$nperblock -lopenblas -g -O3 -o main.out
cp ./main.out $dir_local
cd ..
# copy necessary files to the tmp directory on node
cp -r ./data $dir_local
cp -r ./cp_middle.sh $dir_local
cp -r ./code $dir_local
cp -r ./model $dir_local
cd $dir_local #now we are utilizing the SSD disk locate at compute node
export OPENBLAS_NUM_THREADS=1

# ------------------------
TEMPLATE=model/EPSOC_z_pf.nml        # your source file
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
mpirun -np "$ntasks" ./main.out > main.log 

# copy file and remove tmp directory on node
rm -r namelists
# post process
cd ./code
gfortran outputnew.f90 -cpp -DMPI -fcheck=all -g -o pp.out
cp pp.out ../
cd ..
./pp.out > pp.log
# copy the whole directory back to storage directory
cd /data
cp -a $dir_local $output_dir
rm -rf $dir_local
