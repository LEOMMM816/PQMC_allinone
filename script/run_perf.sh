#!/bin/bash

# 针对 OpenMPI (gfortran通常搭配这个)
RANK=$OMPI_COMM_WORLD_RANK

# 如果你用的是 MPICH 或 Intel MPI，取消下面这行的注释并注释掉上面那行
# RANK=$PMI_RANK

# 只有 Rank 0 运行 perf，其他进程正常运行
if [ "$RANK" -eq "0" ]; then
    echo "Rank $RANK: Running with perf..."
    # -o 指定输出文件名，避免混乱
    perf record -g -o perf.data.rank0 "$@"
else
    # "$@" 代表传递给脚本的所有参数（即你的程序名）
    "$@"
fi