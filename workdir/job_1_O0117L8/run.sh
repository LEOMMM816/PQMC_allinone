
    cd workdir/job_1_O0117L8
    echo "Running in workdir/job_1_O0117L8"
    # 每次运行前清理环境，确保纯净
    
echo "=============================="

    # 环境变量
    source ./param.env
    # 运行公共程序
     mpirun -np 8 /home/zhangxiangyu/project/all_in_one_git/bin/simulation_app > main.log
    # 计算部分完成，开始数据处理
     mpirun -np 1 /home/zhangxiangyu/project/all_in_one_git/bin/postpro_app > pp.log
    