#!/bin/bash

# i need a for loop to sbatch multiple jobs with different parameters
# in each loop, copy model/params.nml. Rename and use sed to modify the copied .nml file
# then sbatch the submit.sh script with the modified .nml file
# jobs_param.tsv is not used here

# include date in the job name to make it unique
JOBNMAE="$(date +%Y%m%d_%H%M%S)_epzOi"

for job in {1..4}

NMLNAME="model/params_job${job}.nml"
cp model/params.nml $NMLNAME
sed -i "s/lat_size = .*/lat_size = $((4 + job*2)),$((4 + job*2))/" $NMLNAME
sed -i "s/beta = .*/beta = $((64 + job*8))d0/" $NMLNAME
#sed -i "scb /filling = .*/filling = $(echo "scale=2; 0.5 + job*0.1" | bc)d0/" $NMLNAME
SBATCH --job-name="${JOBNMAE}_job${job}" \
       --output="logs/${JOBNMAE}_job${job}.out" \
       --error="logs/${JOBNMAE}_job${job}.err" \
       --partition=intel-sc3-32c \
       --qos=normal \
       --nodes=1 \
       --ntasks-per-node=64 \
       --mem-per-cpu=5G \
       script/mpi_32c.sh
done

# sbatch -j epzL10Oi -p intel-sc3-32c -q normal -N 1 --ntasks-per-node=64  ./run_batch.sh -- 1

# sbatch -j epzL10Oi -p intel-sc3-32c -q normal -N 1 --ntasks-per-node=64  ./run_batch.sh -- 2

