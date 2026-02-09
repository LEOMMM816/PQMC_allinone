# %%
# open reading_guide.nml to extract what's in out_core*.dat files
# first read in the basics in reading guide
# define a function that can read namelist format files from fortran
import numpy as np
import os
import re
from pyparsing import line
import matplotlib.pyplot as plt
def read_namelist(filename, namelist_name):
    with open(filename, 'r') as f:
        lines = f.readlines()
        for i, line in enumerate(lines):
            if line.strip().startswith('&' + namelist_name):
                namelist_data = {}
                for j in range(i + 1, len(lines)):
                    if lines[j].strip().startswith('/'):
                        break
                    match = re.match(r'(\w+)\s*=\s*(.*)', lines[j].strip())
                    if match: 
                        key, value = match.groups()
                        namelist_data[key] = eval(value)
                return namelist_data
def read_obs_info(filename):
# read observable info from reading_guide.nml
# there are multiple &obs namelists, read them all and return a dict[obs_name] = {info}
    with open(filename, 'r') as f:
        lines = f.readlines()
        flag = False
        obs_info = {}
        current_obs = {}
        for i, line in enumerate(lines):
            if line.strip().startswith('&obs'):
                flag = True
                current_obs = {}
            if line.strip().startswith('/') and current_obs:
                flag = False
                if 'name' in current_obs:
                    obs_info[current_obs['name']] = current_obs
            if flag:
                match = re.match(r'(\w+)\s*=\s*(.*)', lines[i].strip())
                if match:
                    key, value = match.groups()
                    current_obs[key] = eval(value)
    return obs_info

def recalculate_k_data_function(data, basic_info, obs_info):
    n_bin = basic_info['n_bins']
    sys_size = basic_info['field_width']
    for obs in obs_info.keys():
        obs_dict = obs_info[obs]
        if obs_dict['kind']==1 or obs_dict['name'].endswith('_k'):
            continue
        r_lo = obs_dict['offset_lo']
        r_hi = obs_dict['offset_hi']
        k_lo = obs_info[f'{obs_dict["name"]}_k']['offset_lo']
        k_hi = obs_info[f'{obs_dict["name"]}_k']['offset_hi']
        r_data = data[r_lo-1:r_hi,:]
        
        k_data = transform_to_k_space(r_data,int(np.sqrt(sys_size)))
            # check if k_data has big imaginary part
        if np.any(np.abs(k_data.imag) > 1e-6):
            print(f"Warning: Large imaginary part detected in k_data for bin {bin} of observable {obs_dict['name']}_k")
            # stop the program or handle the issue as needed
            print(np.max(np.abs(k_data.imag)))
            # raise SyntaxError("Large imaginary part detected in k_data")
        data[k_lo-1:k_hi,:] = np.real(k_data)
    return data


def transform_to_k_space(data_matrix, L):
    """
    针对形状为 (L*L, M) 的数据进行批量二维 FFT 变换。
    利用“取实部”技巧代替手动对称化。

    参数:
    data_matrix: 2D numpy array, shape = (L*L, 样本数 M)
                 每一列是一个样本。
    L: int, 晶格线性尺寸 (例如 12)

    返回:
    k_data: 2D numpy array, shape = (L*L, 样本数 M)
            全实数数据。
    """
    # 1. 批量重塑
    # 输入是 (L*L, M)，我们要把它变成 (L, L, M)
    # 使用 order='F' 确保先填充第一维(x)，再填充第二维(y)，最后才是样本维
    batch_grid = data_matrix.reshape((L, L, -1), order='F')
    
    # 【注意】现在的空间维度是轴 0 和 轴 1，样本维度是轴 2
    spatial_axes = (0, 1)

    # 2. 预处理 shift (ifftshift)
    # 将几何中心移到 (0,0)
    batch_shifted = np.fft.ifftshift(batch_grid, axes=spatial_axes)
    
    # 3. 批量 FFT
    # 对前两个维度（空间维度）做变换，并行处理所有样本
    k_batch_complex = np.fft.fft2(batch_shifted, axes=spatial_axes)
    
    # 4. 【关键优化】取实部
    # 等价于在实空间做了 f(r) = (f(r) + f(-r))/2 的对称化操作
    k_batch_real = k_batch_complex.real/(L*L)
    
    # 5. 后处理 shift (fftshift)
    # 将 k=0 移回中心
    k_centered = np.fft.fftshift(k_batch_real, axes=spatial_axes)
    
    # 6. 展平返回
    # 将 (L, L, M) 变回 (L*L, M)
    return k_centered.reshape((-1, data_matrix.shape[1]), order='F')

def get_obs_data(obs,block_id,jn,mean_data,err_data,Lerr_data,basic_data, obs_info_data):
    re = {}
    obs_info= obs_info_data[jn]
    basic_info = basic_data[jn]
    if(not obs in obs_info.keys()):
        print(f"Observation {obs} not found in {jn}")
        raise KeyError
    obs_dict = obs_info[obs]
    lo = obs_dict['offset_lo']
    hi = obs_dict['offset_hi']
    re['val'] = mean_data[jn][lo-1:hi,block_id] # [field_with,nblock]
    re['err'] = err_data[jn][lo-1:hi,block_id]
    re['Lerr'] = Lerr_data[jn][lo-1:hi,block_id]
    return re
        

# %%
from pathlib import Path
import os

if "__file__" in globals():
    # 在正常脚本 / Run Python File / Run Current File in Interactive Window 时
    ROOT = Path(__file__).resolve().parent
else:
    # 在纯 notebook 环境（__file__ 不存在）时
    ROOT = Path.cwd()

os.chdir(ROOT)  # 以后所有相对路径都以这个目录为基准
print(f"Current working directory set to: {ROOT}")
# 下面就可以放心用相对路径了

# %%
work_dir = '../zxy_work_dir_20260123/'
shared_jn = 'O0124'
recalculate_k_data = True
all_data = {}
mean_data = {}
err_data = {}
Lerr_data = {}
bmean_data = {}
basic_data = {}
obs_info_data = {}
jn_list=[]
for file in os.listdir(work_dir):
    #print(f'Processing directory: {file}')
    # jobname要截取file中字母L以后的内容,包含L
    jobname = 'L' + file.split('L', 1)[-1]
    jn_list.append(jobname)
    print(f'Jobname extracted: {jobname}')
    basic_info = read_namelist(f'{work_dir}/{file}/data/data/reading_guide.nml', 'basic')
    # print(basic_info.keys())
    obs_info= read_obs_info(f'{work_dir}/{file}/data/data/reading_guide.nml')
    basic_data[jobname] = basic_info
    obs_info_data[jobname] = obs_info
    #print(obs_info['BFx_BFx_k'])
    n_core = basic_info['n_cores']
    n_block = basic_info['MPI_nblock']
    core_per_block = basic_info['MPI_one_block']
    n_bin = basic_info['n_bins'] 
    n_obs = basic_info['n_obs']
    total_width = basic_info['total_width']

    # extract core data into all_data
    for i_core in range(n_core):
        core_file = f'{work_dir}/{file}/data/data/out_core{i_core}.bin'
        if(not Path(core_file).exists()):
            continue
        data = np.fromfile(core_file, dtype=np.float64) # mean(total_with,nbin)
        data = data.reshape((total_width, n_bin), order='F')
        if(recalculate_k_data):
            data = recalculate_k_data_function(data,basic_info, obs_info)
    # now the data contains mean-value of all bins and all obs in one core
    # put the data into all_data
        b_id = int(i_core/core_per_block)
        core_id_in_block = i_core % core_per_block
        if jobname not in all_data:
            all_data[jobname] = {}
        if b_id not in all_data[jobname]:
            all_data[jobname][b_id] = {}
        all_data[jobname][b_id][core_id_in_block] = data
    
    # calculate the err and Lerr for this job
    mean_data[jobname] = np.zeros((total_width,n_block))
    bmean_data[jobname] =  np.zeros((total_width,n_core,n_block))
    err_data[jobname] = np.zeros((total_width,n_block))
    Lerr_data[jobname] = np.zeros((total_width,n_block))
    # sum and divide for mean
    for b_id in range(n_block):
        for i_core in range(core_per_block):
            for bin in range(n_bin):
                mean_data[jobname][:,b_id] = mean_data[jobname][:,b_id] + all_data[jobname][b_id][i_core][:,bin]
                bmean_data[jobname][:,i_core,b_id] = bmean_data[jobname][:,i_core,b_id] + all_data[jobname][b_id][i_core][:,bin]
    mean_data[jobname] = mean_data[jobname] / (core_per_block * n_bin)
    bmean_data[jobname] = bmean_data[jobname] / n_bin
    # square, sum and divide for err
    for b_id in range(n_block):
        for i_core in range(core_per_block):
            Lerr_data[jobname][:,b_id] += (bmean_data[jobname][:,i_core,b_id] - mean_data[jobname][:,b_id])**2
            for i_bin in range(n_bin):
                err_data[jobname][:,b_id] += (all_data[jobname][b_id][i_core][:,i_bin] - mean_data[jobname][:,b_id])**2
        Lerr_data[jobname][:,b_id] = np.sqrt(Lerr_data[jobname][:,b_id] / core_per_block/(core_per_block-1))
        err_data[jobname][:,b_id] = np.sqrt(err_data[jobname][:,b_id] / (core_per_block * n_bin) /(core_per_block * n_bin-1))
    
# %%
dump_data ={}
for jn in jn_list:
    if('BC' in jn):
        jn_root='BC' + shared_jn
    else:
        jn_root=shared_jn
    jn_temp = jn.split('L')[1]
    if shared_jn in jn_temp:
        # get rid of shared_jn
        jn_temp = jn_temp.replace(shared_jn, '')
    if 'BC' in jn_temp:
        jn_temp = jn_temp.replace('BC', '')
    Lsize = int(''.join(filter(str.isdigit, jn_temp)))
    print(f"Processing jobname: {jn}, Lsize: {Lsize}, jn_root: {jn_root}")
    for b_id in range(basic_data[jn]['MPI_nblock']):
        for obs in obs_info_data[jn].keys():
            # print(f"Processing block: {b_id}, observation: {obs}")
            temp = get_obs_data(obs,b_id,jn,mean_data,err_data,Lerr_data,basic_data, obs_info_data)
            for k,v in temp.items():
                # turn all np object in to list
                if isinstance(v, np.ndarray):
                    temp[k] = v.tolist()
            dump_data.setdefault(jn_root, {}).setdefault('L'+ str(Lsize), {}).setdefault('b'+ str(b_id+1), {}).setdefault(obs, {}).update(temp)
            
print(dump_data[shared_jn]['L12'].keys())
# %%
import json
output_path = ROOT/"extracted_data.json"
with open(output_path, "w") as f:
    json.dump(dump_data, f)
print(f"Data extracted and saved to {output_path}")
# now we can load the data from json for later use

# %%
