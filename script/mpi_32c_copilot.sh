#!/bin/bash
#SBATCH --job-name=${SLURM_JOB_NAME:-pqmc}
#SBATCH --output=/storage/wucongjunLab/zhangxiangyu/tmp/slurm-%j.out
#SBATCH --error=/storage/wucongjunLab/zhangxiangyu/tmp/slurm-%j.err

#set -euo pipefail
#set -x
# ----------------- [编译前检查] -----------------
echo "================= 节点硬件检查 ================="
grep -o 'avx512[^ ]*' /proc/cpuinfo | sort -u | head -n 1
gcc -march=native -Q --help=target | grep march
echo "================= 检查结束 ====================="

mkdir -p /storage/wucongjunLab/zhangxiangyu/tmp
# ---------------------------
# Robust SLURM environment info (use SLURM-provided vars where possible)
submit_dir=${SLURM_SUBMIT_DIR:-$(pwd)}
jobname=${SLURM_JOB_NAME:-job_${SLURM_JOB_ID:-unknown}}
submit_dir=${submit_dir%/}   # strip trailing slash
output_dir="${submit_dir}/results/${jobname}"

# Prefer SLURM_NNODES and SLURM_NTASKS (more standard)
nodes=${SLURM_NNODES:-${SLURM_JOB_NUM_NODES:-1}}
ntasks=${SLURM_NTASKS:-0}

# derive tasks_per_node robustly
if [[ -n "${SLURM_TASKS_PER_NODE:-}" ]]; then
    tasks_per_node=${SLURM_TASKS_PER_NODE%%(*}   # gets number before parentheses like "64(x2)"
elif [[ $nodes -gt 0 && $ntasks -gt 0 ]]; then
    tasks_per_node=$(( ntasks / nodes ))
else
    tasks_per_node=1
fi

# NBLOCK should be exported from the python submitter (or provide a default)
NBLOCK=${NBLOCK:-1}
# fallback ntasks calculation if ntasks not provided by SLURM
if [[ $ntasks -eq 0 ]]; then
    ntasks=$(( nodes * tasks_per_node ))
fi

# node-local temp directory (use jobid if available)
jobid=${SLURM_JOB_ID:-$$}
dir_local="${submit_dir}/${jobname}_${jobid}"
mkdir -p "$dir_local"
chmod 700 "$dir_local"

echo "===== SLURM INFO ====="
echo "submit_dir=$submit_dir"
echo "jobname=$jobname"
echo "jobid=$jobid"
echo "nodes=$nodes"
echo "ntasks=$ntasks"
echo "tasks_per_node=$tasks_per_node"
echo "NBLOCK=$NBLOCK"
echo "output_dir=$output_dir"
echo "dir_local=$dir_local"

ulimit -c unlimited || true
ulimit -a || true

module purge
module load gcc/9.4.0
module load openmpi/4.1.0
module load OpenBLAS/0.3.13

echo "module list (below):"
module list || true

echo "which mpirun: $(which mpirun || true)"
mpirun --version || true
srun --version || true

export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1

# ------------------------

cd "${submit_dir}"
# copy executable and necessary files to node-local dir
cp -r "${submit_dir}/script" "${dir_local}/" || true
cp -r "${submit_dir}/code" "${dir_local}/"
mkdir -p "${dir_local}/model"
cp -r "${submit_dir}/model/${jobname}" "${dir_local}/model/" || true
# switch to node-local dir for runtime
cd "${dir_local}"
mkdir -p ./data/data ./data/out_files
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
# compile the fortran code (compile in the node-local dir to avoid NFS issues)
cd ./code
echo "----- Compiling MPI program -----"
rm -rf ./*.mod
rm -rf ./*.out
mpifort mod_nrtype.f90 input.f90  mod_nrutil.f90 mod_matrixlib.f90 mod_ranstate.f90 mod_lattice.f90 mod_phonon_field.f90 mod_evolution.f90 mod_update.f90 mod_meas.f90 Main_PQMC.f90 \
 -cpp -DMPI -DCPPBLOCK=$NBLOCK -DCPPPERBLOCK=$(( ntasks / NBLOCK )) -march=icelake-server -lopenblas -g -O3 -o main.out
cp ./main.out ../
cd ..
echo "================= 二进制文件检查 ================="
# 2. 确认你的程序确实编译出了 AVX-512 指令 (zmm 寄存器)
# 如果 main.out 里有 zmm，这里会打印出来；如果没有，则是空的
if [ -f "main.out" ]; then
    objdump -d main.out | grep zmm | head -n 5
else
    echo "Error: main.out not found!"
fi
echo "================= 检查结束 ====================="
echo "----- Running MPI program -----"
echo "Executable path: $(pwd)/main.out"

# Recommended: use srun (SLURM-native). This usually avoids mpirun/pmix/ssh integration issues.
# Write both stdout and stderr to main.log
mkdir -p "${dir_local}/run_logs"
mpirun --display-map --map-by ppr:${tasks_per_node}:node -np "$ntasks" --bind-to core --report-bindings --tag-output  ./main.out 2>&1 | tee -a "${dir_local}/run_logs/main_mpirun.log" || {
    echo "mpi failed or returned non-zero. Trying mpirun fallback..."
    # Fallback for mpirun (only use if OpenMPI is SLURM-aware). Try a conservative mapping.
    srun --mpi=pmi2 -n "$ntasks" --ntasks-per-node="$tasks_per_node" ./main.out &>> main.log || {
        echo "Both srun and mpirun failed. Exiting with error."
        exit 1
    }
}

echo "----- post processing -----"
cd ./code
gfortran outputnew.f90 -cpp -DMPI -fcheck=all -g -o pp.out || true
cp pp.out ../ || true
cd ..
./pp.out &> pp.log || true

# ensure output directory exists on the shared filesystem (use absolute path)
mkdir -p "${output_dir}"

if cp -a "${dir_local}/." "${output_dir}/"; then
    echo "cp 成功"
    rm -rf "$dir_local"
else
    echo "cp 失败，尝试用 rsync 作为备用方案..."
    if rsync -a "${dir_local}/" "${output_dir}/"; then
        echo "rsync 成功"
        rm -rf "$dir_local"
    else
        echo "rsync 也失败了，退出脚本"
        exit 1
    fi
fi

echo "Job finished at $(date)"