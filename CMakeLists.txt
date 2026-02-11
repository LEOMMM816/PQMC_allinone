cmake_minimum_required(VERSION 3.10)
project(HPC_Project LANGUAGES Fortran C CXX)

# =========================================================
#  自动生成的 CMakeLists.txt
#  编译器模式: GNU
#  生成时间: Wed Feb 11 21:03:37 CST 2026
# =========================================================

# 设置 C++ 标准 (如果用到)
set(CMAKE_CXX_STANDARD 14)

# ---------------------------------------------------------
# 1. 编译器 Flag 配置
# ---------------------------------------------------------

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    message(STATUS ">>> 配置 GNU (gfortran) 环境")
    
    # Release 模式 Flags
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -g -march=native -ffree-line-length-none -fbacktrace -cpp -DMPI")
    
    # Debug 模式 Flags
    set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -Wall -fcheck=all -fbacktrace -cpp -DMPI -ffpe-trap=invalid,zero,overflow")
    
    # set(CMAKE_Fortran_FLAGS "-O3 -g -march=native -ffree-line-length-none -fbacktrace -cpp -DMPI")
    # set(CMAKE_Fortran_FLAGS "-O0 -g -Wall -fcheck=all -fbacktrace -cpp -DMPI -ffpe-trap=invalid,zero,overflow")
    set(CMAKE_Fortran_FLAGS "-O3 -g -march=native -ffree-line-length-none -fbacktrace -cpp -DMPI")
            
endif()
# ---------------------------------------------------------
# 2. MPI 配置
# ---------------------------------------------------------

# --- MPI 配置 ---
find_package(MPI REQUIRED)
# 兼容旧版 CMake 的宏定义写法
add_definitions(-DMPI)

# 如果 Fortran 编译器就是 MPI 包装器 (如 mpiifort)，
# 上面的 find_package 主要是为了保险和查找头文件。
include_directories(${MPI_Fortran_INCLUDE_PATH})

# --- 自动归纳文件位置 ---
# 1. 可执行文件 -> 项目根目录/bin
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

# 2. Fortran .mod 文件 -> 项目根目录/modules
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_SOURCE_DIR}/modules)

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
        
if(MY_MATH_LIBS)
    target_link_libraries(core_mods PUBLIC ${MY_MATH_LIBS})
endif()
if(MPI_FOUND)
    target_link_libraries(core_mods PUBLIC MPI::MPI_Fortran)
endif()
target_link_libraries(simulation_app PRIVATE core_mods)
target_link_libraries(postpro_app PRIVATE core_mods)

