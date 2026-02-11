import sys
import shutil
import re
import os
from pathlib import Path
from datetime import date
import subprocess  
import argparse
# ----------------- 对象 ----------------- #
# 创建一个对象来存储任务参数，task对象包含以下属性：
# - id: 任务的唯一标识符
# - template: 使用的模板文件名
# - sbatch_param: 一个字典，存储所有以 "sbatch-" 开头的参数
# - nml_param: 一个字典，存储所有非 "sbatch-" 和 "env-" 开头的参数(除了template)
# - env_param: 一个字典，存储所有以 "env-" 开头的参数
class Task:
    def __init__(self, id, template):
        self.id = id
        self.template = template
        self.sbatch_param = {}
        self.nml_param = {}
        self.env_param = {}



# ----------------- 工具函数 ----------------- #
def parse_value(raw):
    """
    将单元格字符串解析成描述字典：
      - scalar: 普通数字/字符串
      - tuple: "(8,8,8)"  -> inner = "8,8,8"
      - list:  "[1,2,3]"
      - range: "[1:0.1:1.7]"
    """
    raw = raw.strip()
    if raw.startswith("(") and raw.endswith(")"):
        inner = raw[1:-1].strip()
        return {"kind": "tuple", "inner": inner}

    if raw.startswith("[") and raw.endswith("]"):
        inner = raw[1:-1].strip()
        # 形如 start:step:stop 的 range
        if inner.count(":") == 2 and "," not in inner:
            start_s, step_s, stop_s = [s.strip() for s in inner.split(":")]
            return {
                "kind": "range",
                "start": float(start_s),
                "step": float(step_s),
                "stop": float(stop_s),
            }
        else:
            # 普通 list: 按逗号切
            parts = [p.strip() for p in inner.split(",") if p.strip()]
            return {"kind": "list", "values": parts}

    # 默认：标量
    return {"kind": "scalar", "value": raw}


def format_float(x):
    # 用通用格式，避免 1.7000000000002 这种
    return f"{x:g}"


def expand_for_blocks(desc, nblock, varname):
    """
    根据 nblock 展开到每个 block 的取值列表（长度 = nblock）。
    desc 为 parse_value 的结果。
    """
    kind = desc["kind"]

    if kind == "scalar":
        return [desc["value"]] * nblock

    if kind == "tuple":
        # 去掉括号后的内容整体写入
        return [desc["inner"]] * nblock

    if kind == "list":
        values = desc["values"]
        if len(values) != nblock:
            raise ValueError(
                f"变量 '{varname}' 的 list 元素数为 {len(values)}，"
                f"但 nblock = {nblock}，不一致，程序终止。"
            )
        return values

    if kind == "range":
        start = desc["start"]
        step = desc["step"]
        stop = desc["stop"]
        if step == 0:
            raise ValueError(f"变量 '{varname}' 的 range 步长为 0，非法。")

        # 先计算理论元素个数
        expected_n = int(round((stop - start) / step)) + 1
        values = []
        for i in range(expected_n):
            val = start + i * step
            # 防止数值误差轻微越界
            if step > 0 and val > stop + 1e-9:
                break
            if step < 0 and val < stop - 1e-9:
                break
            values.append(format_float(val))

        if len(values) != nblock:
            raise ValueError(
                f"变量 '{varname}' 的 range [{start}:{step}:{stop}] "
                f"生成了 {len(values)} 个元素，但 nblock = {nblock}，不一致，程序终止。"
            )
        return values

    raise ValueError(f"未知 kind: {kind}")


def replace_in_nml(text, varname, new_value):
    """
    在 namelist 文本中，把一行形如
      varname = 原值
    改成
      varname = new_value
    如果找不到该变量，则在末尾追加一行。
    """
    pattern = re.compile(rf"^(\s*{re.escape(varname)}\s*=\s*).*$", re.MULTILINE)
    if pattern.search(text):
        # 用函数替换，避免 \1{new_value} 变成 \18 之类的问题
        def _repl(m):
            return m.group(1) + str(new_value)
        new_text = pattern.sub(_repl, text, count=1)
        return new_text
    else:
        # 没找到变量，就追加到文件末尾
        print(f"警告: 在 nml 文件中找不到变量 '{varname}'，将在末尾追加该变量。", file=sys.stderr)
        append_line = f"\n{varname} = {new_value}\n"
        return text + append_line
    
def convert_tsv_lines_to_dicts(lines):

    """
    将 TSV 文件的多行内容转换为字典列表。
    第一行是表头，后续行是数据。
    忽略空行和以 ! 开头的注释行。
    """
    headers = [h.strip() for h in lines[0].strip().split()]
    # 检查必须列
    if "template" not in headers:
        print("表头中找不到 'template' 列。", file=sys.stderr)
        sys.exit(1)
    if "sbatch-job-name" not in headers:
        print("表头中没有 'sbatch-job-name' 列。", file=sys.stderr)
        sys.exit(1)
    if "env-nblock" not in headers:
        print("表头中没有 'env-nblock' 列。", file=sys.stderr)
        sys.exit(1)
    tasks = []
    for line in lines[1:]:
        if not line.strip() or line.lstrip().startswith("!"):
            continue  # 忽略空行和注释行
        parts = [p.strip() for p in line.strip().split()]
        if len(parts) != len(headers):
            raise ValueError("TSV 行的列数与表头不匹配。")
        temp_dict = dict(zip(headers, parts))
        #开始读取task内容 "***_param"'s header takes the form of "***-headername"
        
        task = Task(id=len(tasks) + 1, template=temp_dict["template"])
        task.sbatch_param = {}
        task.nml_param = {}
        task.env_param = {}
        for key, value in temp_dict.items():
            if key == "template":
                continue  # 已经处理过 template 了
            if key.startswith("sbatch-"):
                param_name = key[len("sbatch-"):]
                task.sbatch_param[param_name] = value
            elif key.startswith("env-"):
                param_name = key[len("env-"):]
                task.env_param[param_name] = value
            else:
                task.nml_param[key] = value
        tasks.append(task)

    return  tasks




#----------------- 生成 CMakeLists.txt ----------------- #
def generate_cmake(project_name, compiler, use_mpi,platform,debug):
    # ---------------------------------------------------------
    # 1. 定义不同编译器的 Flag 配置 (核心逻辑)
    # ---------------------------------------------------------
    
    # === Intel 编译器配置 (针对集群环境优化) ===
    # -mkl: 自动链接数学库 (兼容旧版写法)
    # -heap-arrays 1024: 防止栈溢出 (Segfault 杀手)
    # -assume byterecl: I/O 兼容性
    # -traceback: 报错时打印行号
    intel_flags_release = ("-O3 -xICELAKE-SERVER -qopt-zmm-usage=high "
    "-mkl -traceback -heap-arrays 1024 -fpp -DMPI -assume byterecl")
    intel_flags_debug   = "-O0 -g -traceback -check all -warn all -fpe0 -mkl -heap-arrays 1024 -fpp -DMPI"

    # === GNU (gfortran) 配置 ===
    # -ffree-line-length-none: 防止代码行过长报错
    # -fbacktrace: 报错打印堆栈
    if(platform=="cluster" or platform=="SSD"):
        # 针对集群的优化
        gnu_flags_release = "-O3 -Ofast -g -march=icelake-server -ffree-line-length-none -fbacktrace -cpp -DMPI"
    else:
        # 本地编译
        gnu_flags_release = "-O3 -g -march=native -ffree-line-length-none -fbacktrace -cpp -DMPI"
    
    gnu_flags_debug   = "-O0 -g -Wall -fcheck=all -fbacktrace -cpp -DMPI -ffpe-trap=invalid,zero,overflow"

    # ---------------------------------------------------------
    # 2. 构建 CMakeLists.txt 内容
    # ---------------------------------------------------------
    content = f"""cmake_minimum_required(VERSION 3.10)
project({project_name} LANGUAGES Fortran C CXX)

# =========================================================
#  自动生成的 CMakeLists.txt
#  编译器模式: {compiler.upper()}
#  生成时间: {os.popen('date').read().strip()}
# =========================================================

# 设置 C++ 标准 (如果用到)
set(CMAKE_CXX_STANDARD 14)

# ---------------------------------------------------------
# 1. 编译器 Flag 配置
# ---------------------------------------------------------
"""

    if compiler == 'intel':
        content += f"""
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    message(STATUS ">>> 配置 Intel 编译器环境")
    
    # Release 模式 Flags
    set(CMAKE_Fortran_FLAGS_RELEASE "{intel_flags_release}")
    
    # Debug 模式 Flags (cmake -DCMAKE_BUILD_TYPE=Debug ..)
    set(CMAKE_Fortran_FLAGS_DEBUG "{intel_flags_debug}")
    
    # 默认通用 Flags
    """
        if(debug):
            content += f"""set(CMAKE_Fortran_FLAGS "{intel_flags_debug}")
            """
        else:
            content += f"""set(CMAKE_Fortran_FLAGS "{intel_flags_release}  -qopt-report=5 -qopt-report-phase=vec")
            """
        content += f"""
else()
    message(WARNING "你选择了生成 Intel 配置，但 CMake 检测到的编译器不是 Intel！")
endif()
"""
    elif compiler == 'gnu':
        content += f"""
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    message(STATUS ">>> 配置 GNU (gfortran) 环境")
    
    # Release 模式 Flags
    set(CMAKE_Fortran_FLAGS_RELEASE "{gnu_flags_release}")
    
    # Debug 模式 Flags
    set(CMAKE_Fortran_FLAGS_DEBUG "{gnu_flags_debug}")
    
    # set(CMAKE_Fortran_FLAGS "{gnu_flags_release}")
    # set(CMAKE_Fortran_FLAGS "{gnu_flags_debug}")
    """
        if(debug):
            content += f"""set(CMAKE_Fortran_FLAGS "{gnu_flags_debug}")
            """
        else:
            content += f"""set(CMAKE_Fortran_FLAGS "{gnu_flags_release}")
            """
        content += f"""
endif()
# ---------------------------------------------------------
# 2. MPI 配置
# ---------------------------------------------------------
"""
    if use_mpi:
        content += """
# --- MPI 配置 ---
find_package(MPI REQUIRED)
# 兼容旧版 CMake 的宏定义写法
add_definitions(-DMPI)

# 如果 Fortran 编译器就是 MPI 包装器 (如 mpiifort)，
# 上面的 find_package 主要是为了保险和查找头文件。
include_directories(${MPI_Fortran_INCLUDE_PATH})
"""
    else:
        content += """
# --- 串行模式 (无 MPI) ---
message(STATUS ">>> 禁用 MPI")
"""
    # 输出文件目的地
    content += f"""
# --- 自动归纳文件位置 ---
# 1. 可执行文件 -> 项目根目录/bin
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${{CMAKE_SOURCE_DIR}}/bin)

# 2. Fortran .mod 文件 -> 项目根目录/modules
set(CMAKE_Fortran_MODULE_DIRECTORY ${{CMAKE_SOURCE_DIR}}/modules)
"""
    content += f"""
# ---------------------------------------------------------
# 3. 构建目标
# ---------------------------------------------------------

add_library(core_mods
    src/mod_nrtype.f90
    src/mod_nrutil.f90
    src/mod_matrixlib.f90
    src/mod_lattice.f90
    src/input.f90
    src/input.f90
    src/mod_nrtype.f90
    src/mod_nrutil.f90
    src/mod_ranstate.f90
    src/mod_matrixlib.f90
    src/mod_lattice.f90
    src/mod_phonon_field.f90
    src/mod_evolution.f90
    src/mod_update.f90
    src/mod_meas.f90
    # ... 其他文件 ...
)
# 主文件
add_executable(simulation_app src/Main_PQMC.f90)
# 后处理文件
add_executable(postpro_app src/outputnew.f90)

# ---------------------------------------------------------
# 4. 库链接 (Link Libraries)
# ---------------------------------------------------------
"""
# GNU 编译器通常需要显式链接 BLAS/LAPACK
    if compiler == 'intel':
        content += """
message(STATUS ">>> Intel 编译器使用 MKL 自动链接数学库")
set(MY_MATH_LIBS "")
""" 

    elif compiler == 'gnu':
        content += """
message(STATUS ">>> 尝试链接 OpenBLAS/LAPACK (GNU 环境)")
 # 查找库
    find_package(BLAS REQUIRED)
    find_package(LAPACK REQUIRED)

    # ---> 输出 BLAS 检查信息 <---
    message(STATUS "    >>> [BLAS 库检查]")
    if(BLAS_FOUND)
        message(STATUS "        状态: 已找到")
        message(STATUS "        库文件位置: ${BLAS_LIBRARIES}")
    endif()

    # ---> 输出 LAPACK 检查信息 <---
    message(STATUS "    >>> [LAPACK 库检查]")
    if(LAPACK_FOUND)
        message(STATUS "        状态: 已找到")
        message(STATUS "        库文件位置: ${LAPACK_LIBRARIES}")
    endif()

    set(MY_MATH_LIBS ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
        """


    # 可执行文件和公共模块、mpi模块链接
    content += """
if(MY_MATH_LIBS)
    target_link_libraries(core_mods PUBLIC ${MY_MATH_LIBS})
endif()
if(MPI_FOUND)
    target_link_libraries(core_mods PUBLIC MPI::MPI_Fortran)
endif()
target_link_libraries(simulation_app PRIVATE core_mods)
target_link_libraries(postpro_app PRIVATE core_mods)

""" 

    # 写入文件
    with open("CMakeLists.txt", "w") as f:
        f.write(content)
    
    print(f"✅ 成功生成 CMakeLists.txt (模式: {compiler}, MPI: {use_mpi})")
    # print(f"👉 下一步: mkdir build && cd build && cmake .. && make -j")

# ----------------- 生成提交脚本 ----------------- #
def generate_script(sb_header, job_dir, exe_path, postpro_path, total_tasks,module_str,platform,compiler):
    content = ""
    if(platform=="SSD" or platform=="cluster"):
        content +=f"""#!/bin/bash
{sb_header}
#SBATCH -c 1                     # 每个任务使用 1 个 CPU 核心
#SBATCH --output={job_dir}/slurm-%j.out  
#SBATCH --error={job_dir}/slurm-%j.err  
"""

    # 进入工作目录
    content += f"""
    cd {job_dir}
    echo "Running in {job_dir}"
    # 每次运行前清理环境，确保纯净
    
echo "=============================="
"""
    if(platform=="SSD" or platform=="cluster"):
    # module 加载
        content += f"""
    {module_str}
    echo "=== Slurm Allocation Check ==="
    echo "Total MPI Tasks (核心总数): $SLURM_NTASKS"
    echo "Nodes (节点数): $SLURM_NNODES"
    echo "Node List (节点列表): $SLURM_NODELIST"
    """
        if(compiler=="intel"):
            content += f"""
    # 【新增】强制 Intel Fortran 实时输出，不缓冲
    export FORT_BUFFERED=no
    # 【新增】强制 Intel MPI 实时输出
    export I_MPI_JOB_OUTPUT_BUFFERING=0
    #export I_MPI_HYDRA_BOOTSTRAP=slurm
    # 关闭 Intel MPI 的自作主张
    export I_MPI_PIN=disable
    export I_MPI_PIN_DOMAIN=auto
    # 确保网络协议正确 (回顾之前的报错)
    export I_MPI_FABRICS=shm:tcp
    echo "=== Actual Running Check ==="
    srun --cpu-bind=none hostname | sort | uniq -c

    # 环境变量
    source ./param.env
    # 运行公共程序
    # 此时程序会读取当前目录下的./model/*.nml 文件
    # 并在 ./data/ 目录下读写数据
    # mpirun -np {total_tasks} {exe_path}
    srun --mpi=pmi2 --cpu-bind=none {exe_path} > main.log

    # 计算部分完成，开始数据处理
    # mpirun -np 1 {postpro_path} > pp.log
    srun --mpi=pmi2 --cpu-bind=none {postpro_path} > pp.log
    """
        else: # gnu on cluster
            content += f"""
    # 环境变量
    source ./param.env
    # 运行公共程序
    mpirun -np {total_tasks} {exe_path} > main.log
    # 计算部分完成，开始数据处理
    mpirun -np 1 {postpro_path} > pp.log
    """
    else: # local
        content += f"""
    # 环境变量
    source ./param.env
    # 运行公共程序
     mpirun -np {total_tasks} {exe_path} > main.log
    # 计算部分完成，开始数据处理
     mpirun -np 1 {postpro_path} > pp.log
    """
    return content
# ----------------- 主逻辑 ----------------- #

def submit_main(compiler, platform, jobfile):
    # 读取 参数列表
    # tsv_addr = './input/' + jobfile # 参数文件路径 
    tsv_addr = Path(jobfile)
    if not tsv_addr.exists():
        print(f"找不到参数文件: {tsv_addr}", file=sys.stderr)
        sys.exit(1)
    #删除整个 work_dir 目录
    if(platform!="local" and platform!="SSD" and platform!="cluster"):
        print("第二个参数应为 local 或 SSD 或 cluster", file=sys.stderr)
        sys.exit(1)
    today_str = date.today().strftime("%Y%m%d")
# ==================================================================
# 0.Cmake
# ==================================================================
    print("正在编译程序...")
    subprocess.run("rm -rf ./build", shell=True)
    subprocess.run("rm -rf ./bin", shell=True)
    subprocess.run("mkdir -p bin", shell=True)
    subprocess.run("mkdir -p build", shell=True)
    module_list = []
    if(platform=="cluster" or platform=="SSD"):
        # 集群编译
        # 把所有module命令整合成一个字符串传给 subprocess.run
        if(compiler=="intel"):
            module_list =  ["module purge",
                            "module load intel/intel2019",
                            "module load intelmpi/2019",
                            "module load cmake/3.19.6",
                            "module load intelmkl/2019"]
            FC_STR = "mpiifort"
        elif(compiler=="gnu"):
            module_list =  ["module purge",
                            "module load gcc/9.4.0",
                            "module load openmpi/4.1.0",
                            "module load cmake/3.19.6",
                            "module load OpenBLAS/0.3.13"]
            FC_STR = "mpifort"
        module_str = "\n".join(module_list)
        build_script_content = f"""#!/bin/bash
set -e  # 遇到错误立即停止

# 加载环境
{module_str}
# 打印一下看看
echo "当前加载的模块:"
module list

# 编译
mkdir -p build
cd build
# 加上 FC 确保万无一失
FC={FC_STR} cmake ..
make
"""
        with open("run_compile.sh", "w") as f:
            f.write(build_script_content)
        print("开始运行编译脚本...")
        subprocess.run("bash run_compile.sh", shell=True)
    elif(platform=="local"):      
        subprocess.run("cd build && cmake .. && make && cd ..", shell=True)
    print("编译完成。")
# ==================================================================
# submit
# ==================================================================
# 1. 预先检查：确保程序已经编译好了
    exe_path = os.path.abspath("./bin/simulation_app")
    if not os.path.exists(exe_path):
        print("Error: 请先运行 cmake && make 编译程序！")
        exit(1)
    postpro_path = os.path.abspath("./bin/postpro_app")
    
    if(platform=="SSD"):
        work_dir = Path(f"/data/zxy_workdir_{today_str}")
    elif platform=="local":
        work_dir = Path("./workdir")
    elif platform=="cluster":
        work_dir = Path(f"./zxy_work_dir_{today_str}")
    if work_dir.exists():
        shutil.rmtree(work_dir)
    # 创建 空白work_dir 目录
    work_dir.mkdir(parents=True, exist_ok=True)
    print(f"工作目录: {work_dir}")

    # 2. 为每个任务创建独立的文件夹
    # 假设 task 有一个唯一的 ID，或者用行号
    # 读取 jobs.tsv，忽略空行和以 ! 开头的行
    with tsv_addr.open("r") as f:
        lines = [ln for ln in f if ln.strip() and not ln.lstrip().startswith("!")]
        if not lines:
            print("jobs.tsv 为空。", file=sys.stderr)
            sys.exit(1)
        tasks = convert_tsv_lines_to_dicts(lines) 


    for task in tasks:
        job_dir = os.path.join(work_dir, f"job_{task.id}_{task.sbatch_param.get('job-name','unnamed')}")
        os.makedirs(job_dir, exist_ok=True)

        # 创建 data 目录 (如果程序需要)
        os.makedirs(os.path.join(job_dir, "data","data"), exist_ok=True)
        os.makedirs(os.path.join(job_dir, "data","out_files"), exist_ok=True)
        os.makedirs(os.path.join(job_dir, "data","temp"), exist_ok=True)
        # 创建 model 目录
        job_model_path = os.path.join(job_dir, "model")
        os.makedirs(job_model_path, exist_ok=True)

# 3. 生成该任务特定的nml文件:先把template复制到job_dir/model，并命名为{template}_b{block_id}.nml,一共复制 nblock 份
        template_path = os.path.join("input", "model", task.template + '.nml')
        if not os.path.exists(template_path):
            print(f"Error: 任务 {task.id} 的模板文件不存在: {template_path}", file=sys.stderr)
            continue
        nml_files = []
        for block_id in range(1, int(task.env_param.get("nblock", "1")) + 1):
            nml_filename = f"{os.path.splitext(task.template)[0]}_b{block_id}.nml"
            nml_path = os.path.join(job_model_path, nml_filename)
            shutil.copy2(template_path, nml_path)
            nml_files.append(nml_path)
        # 修改 nml 文件中的参数
        # 解析各个nml变量，并展开到每个 block
        per_var_values = {}  # var -> [v_block1, v_block2, ...]
        for key,var in task.nml_param.items():
            desc = parse_value(var) 
            values = expand_for_blocks(desc, int(task.env_param["nblock"]), var)
            per_var_values[key] = values
        # 为该任务的每个 block 的 nml 文件做替换
        for ib, nml_path in enumerate(nml_files):
            text = Path(nml_path).read_text()
            for var, values in per_var_values.items():
                new_val = values[ib]
                text = replace_in_nml(text, var, new_val)
            Path(nml_path).write_text(text)

        # 设置环境变量文件 param.env
        env_lines = []
        for env_var, env_value in task.env_param.items():
            env_lines.append(f'export {env_var}={env_value}\n')
        env_lines.append(f'export OMP_NUM_THREADS=1\n')  # 固定设置 OMP_NUM_THREADS=1
        env_lines.append(f'export MKL_NUM_THREADS=1\n')
        env_lines.append(f'export OPENBLAS_NUM_THREADS=1\n')
        env_lines.append(f'export I_MPI_HYDRA_BOOTSTRAP=slurm\n')
        #env_lines.append(f'export I_MPI_DEBUG=5\n')  # 开启调试模式，看看卡在哪
        #env_lines.append(f'export I_MPI_HYDRA_BOOTSTRAP=fork\n')
        #env_lines.append(f'export I_MPI_FABRICS=shm:tcp\n')
        env_sh_path = os.path.join(job_dir, "param.env")
        with open(env_sh_path, "w") as f:
            f.writelines(env_lines)

        # 4. 生成提交脚本 run.sh
        # 关键点：脚本里要引用那个公共的 exe_path

        # 设置sbatch变量
        sb_lines = []
        for key, val in task.sbatch_param.items():
            if val is None:
                sb_lines.append(f"#SBATCH --{key}")
            else:
                sb_lines.append(f"#SBATCH --{key}={val}")
        sb_header = "\n".join(sb_lines)
        # 转换module_list为字符串
        module_str = "\n".join(module_list)

        # 提前计算total_tasks
        total_tasks = int(task.sbatch_param.get('nodes', '1')) * int(task.sbatch_param.get('ntasks-per-node', '1'))
        print(f"任务 {task.id} 总任务数: {total_tasks}")
        script_content = generate_script(
            sb_header, job_dir, exe_path, postpro_path, total_tasks,module_str,platform,compiler
        )

        # 写入脚本
        script_path = os.path.join(job_dir, "run.sh")
        with open(script_path, "w") as f:
            f.write(script_content)

        # 5. 提交
        # 因为我们在脚本里写了绝对路径的 cd，所以在哪里提交都行
        if(platform=="local"):
            # 本地运行
            subprocess.run(["bash", script_path])
            print(f"Locally ran job {task.id}")
        elif(platform=="cluster"):
            # 提交到调度系统
            subprocess.run(["sbatch", script_path])
            print(f"Submitted job {task.id}")
        elif(platform=="SSD"):
            # 提交到调度系统，使用SSD存储
            subprocess.run(["sbatch", script_path])
            print(f"Submitted job {task.id} with SSD")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="自动生成适用于 HPC 的 CMakeLists.txt")
    
    parser.add_argument('--compiler', choices=['intel', 'gnu'], default='intel', 
                        help='选择编译器: intel (ifort/ifx) 或 gnu (gfortran)')
    
    parser.add_argument('--mpi', action='store_true', default=True,
                        help='是否启用 MPI (默认启用, 使用 --no-mpi 关闭)')
    parser.add_argument('--no-mpi', action='store_false', dest='mpi',
                        help='关闭 MPI')
    
    parser.add_argument('--name', default='simulation_app',
                        help='生成的可执行文件名 (默认为 simulation_app)')
    
    parser.add_argument('--platform',default='local',choices=['local','cluster','SSD'],
                        help='编译平台选择: local (本地), cluster (集群), SSD (高速存储节点)')

    parser.add_argument('--project', default='HPC_Project',
                        help='CMake 项目名称')
    parser.add_argument('--jobfile', required=True,help='任务参数文件 jobs.tsv 的路径')

    parser.add_argument('--debug', action='store_true', default=False,
                        help='是否启用调试模式 (默认关闭)')

    args = parser.parse_args()
    
    if(args.debug):
        print("⚠️  调试模式已启用！编译器 Flags 将切换到 Debug 模式，且脚本中会设置更多调试环境变量。")
    generate_cmake(args.project, args.compiler, args.mpi, args.platform, args.debug)
    submit_main(args.compiler, args.platform, args.jobfile)
