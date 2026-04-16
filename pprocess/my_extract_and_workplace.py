

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
import re
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
# here I from scratch write an independent version to visionalize data
jobname_list = ["O1644","BCO1644"]
preaddr = "../zxy_work_dir_20260323/"
n_block = 8
# omega in 1.6,1.7...2.3, kept to 2 significant digits
omega = [round(1.6 + 0.4*i,3) for i in range(n_block)]
lat_size = [6, 8, 10,12,14]
data = {}
for jn in jobname_list:
    data[jn] = data.get(jn, {})
    for size in lat_size:
        file_name = f"L{size}{jn}"
        size_name = f"L{size}"
        data[jn][size_name] = data[jn].get(size_name, {})
        # Search under ROOT/preaddr for a folder whose name contains file_name.
        search_root = (ROOT / preaddr).resolve()
        matched_folders = sorted(
            [p for p in search_root.iterdir() if p.is_dir() and file_name in p.name],
            key=lambda p: p.name,
        )
        if not matched_folders:
            print(f"skipped {search_root} as no folder matches '*{file_name}*'")
            continue
        if len(matched_folders) > 1:
            print(f"Multiple folders matched '*{file_name}*', using {matched_folders[0].name}")

        base_path = matched_folders[0] / "data" / "out_files"


        if not base_path.exists():
            print(f"skipped {base_path} as it does not exist")
            continue
        print(f"Loading data from {file_name} ...")
        for block in range(1, n_block+1):
            block_name = f"b{block}"
            data[jn][size_name][block_name] = {}
            # find the file that starts with "b{block}" and ends with ".dat"
            files = list(base_path.glob(f"{block_name}*.out"))
            if not files:
                raise FileNotFoundError(f"No file matching {block_name}*.out in {base_path}")
            file_path = files[0]
            with open(file_path, "r") as f:
                lines = f.readlines()
            obs_begin = False
            for line in lines:
                line = line.strip()
                # if line does not start with "#--observables--#", skip it because observables start from there
                if line.startswith("#--observables--#"):
                    obs_begin = True
                if not obs_begin:
                    continue                
                if line.startswith("#"):
                    continue
                if line == "":
                    continue
                parts = re.split(r'\s+', line)
                name = parts[0]
                data[jn][size_name][block_name][name] = data[jn][size_name][block_name].get(name, {"val": [], "err": [], "Lerr": []})
                try: 
                    val = float(parts[2])
                    err = float(parts[3])
                    Lerr = float(parts[4])
                    data[jn][size_name][block_name][name]["val"].append(val)
                    data[jn][size_name][block_name][name]["err"].append(err)
                    data[jn][size_name][block_name][name]["Lerr"].append(Lerr)
                except ValueError:
                    if name != "Name":
                        print(f"Could not find value for {name} in {jn} {size_name} {block_name}, skipping...")
                    continue

# now we have data, we need to output it to a portable format like json or pickle for later use
import json
output_path = ROOT/"extracted_data.json"
#with open(output_path, "w") as f:
#    json.dump(data, f)
#print(f"Data extracted and saved to {output_path}")
# now we can load the data from json for later use
#with open(output_path, "r") as f:
    #data = json.load(f)
#print(f"Data loaded from {output_path}")


# %% 
# average the data from with BC and without BC
jobname_list_avg = ["O1644"]
for jn in jobname_list_avg:
    jn1 = jn  # without BC
    jn2 = "BC" + jn  # with BC
    jn_target = f"avg_{jn}"
    data[jn_target] = {}
    for size in lat_size:
        size_name = f"L{size}"
        data[jn_target][size_name] = {}
        for block in range(1, n_block+1):
            block_name = f"b{block}"
            data[jn_target][size_name][block_name] = {}
            keys = set(data[jn1][size_name][block_name].keys()).intersection(set(data[jn2][size_name][block_name].keys()))
            for key in keys:
                val1 = np.array(data[jn1][size_name][block_name][key]["val"])
                err1 = np.array(data[jn1][size_name][block_name][key]["err"])
                Lerr1 = np.array(data[jn1][size_name][block_name][key]["Lerr"])
                val2 = np.array(data[jn2][size_name][block_name][key]["val"])
                err2 = np.array(data[jn2][size_name][block_name][key]["err"])
                Lerr2 = np.array(data[jn2][size_name][block_name][key]["Lerr"])
                avg_val = (val1 + val2) / 2
                avg_err = np.sqrt(err1**2 + err2**2) / 2
                avg_Lerr = np.sqrt(Lerr1**2 + Lerr2**2) / 2
                data[jn_target][size_name][block_name][key] = {"val": avg_val.tolist(), "err": avg_err.tolist(), "Lerr": avg_Lerr.tolist()}
# %%
#plot some data to compare different BC
example_name = "den_den"
jn1 = jobname_list_avg[0]  # without BC
jn2 = "BC" + jn1  # with BC
for size in lat_size:
    size_name = f"L{size}"
    block_name = "b8"
    val1 = np.array(data[jn1][size_name][block_name][example_name]["val"])
    err1 = np.array(data[jn1][size_name][block_name][example_name]["err"])
    Lerr1 = np.array(data[jn1][size_name][block_name][example_name]["Lerr"])
    val2 = np.array(data[jn2][size_name][block_name][example_name]["val"])
    err2 = np.array(data[jn2][size_name][block_name][example_name]["err"])
    Lerr2 = np.array(data[jn2][size_name][block_name][example_name]["Lerr"])
    plt.errorbar(size,val1[0], yerr=Lerr1[0], marker='o', linestyle='-', color='blue', label=f'{jn1} {size_name} {block_name} {example_name}')
    plt.errorbar(size,val2[0], yerr=Lerr2[0], marker='o', linestyle='-', color='red', label=f'{jn2} {size_name} {block_name} {example_name}')
    



# %%
# ==========================================
# Data Preprocessing: Normalization, Background, and C(0,0) Subtraction
# ==========================================
import numpy as np

print("Starting data preprocessing...")

for jn in data.keys():
    for size_id in data[jn].keys():
        # Extract integer size L from size_id (e.g., "L14" -> 14)
        L = int(size_id[1:])
        L_square = L**2
        
        for block_id in data[jn][size_id].keys():
            obs_dict = data[jn][size_id][block_id]
            
            # ---------------------------------------------------
            # Step 1: Normalize spin correlations by L^2
            # ---------------------------------------------------
            spin_obs = ["Sx_Sx", "Sy_Sy", "Sx_Sx_k", "Sy_Sy_k"]
            for obs in spin_obs:
                if obs in obs_dict:
                    # Divide values and errors by L^2
                    obs_dict[obs]["val"] = [v / L_square for v in obs_dict[obs]["val"]]
                    obs_dict[obs]["err"] = [e / L_square for e in obs_dict[obs]["err"]]
                    obs_dict[obs]["Lerr"] = [e / L_square for e in obs_dict[obs]["Lerr"]]

        
print("Data preprocessing completed successfully!")
# %%
# plot certain observable vs real(momentum) space for a given lattice size and block(omega)
jn = jobname_list_avg[0]
jn = f"BC{jn}"
size_i = 14
size_id = f"L{size_i}"
block_id = "b8"
xlabel = "x"
#ylabel = "BFxy_BFxy_k"
#ylabel = "SC_SC_k"
ylabel = "Sx_Sx_k"
title = f"{jn}'s {size_id} {ylabel}"
x = range(size_i**2)  # example x values
y = data[jn][size_id][block_id][ylabel]["val"]
#y[int(size_i**2/2 + size_i/2)] = 0  # set the middle point to 0 for better visualization
yerr = data[jn][size_id][block_id][ylabel]["err"]
Lerr = data[jn][size_id][block_id][ylabel]["Lerr"]
plt.errorbar(x, y, yerr=Lerr, marker='o', linestyle='-', label=f'{size_id} {block_id} {ylabel}')
plt.xlabel(xlabel)
plt.ylabel(ylabel)
plt.title(title)
plt.legend()
plt.show()


# %%
# plot observable vs different blocks(omega) for lattice sizes
jn = jobname_list_avg[0]
block_id = [f"b{i}" for i in range(1,9)]
ylabel = ["Sz_Sz_k", "Sx_Sx_k"]
all_sizes = lat_size[1:]  # skip the smallest size for better visualization
all_omegas = [omega[i-1] for i in range(1,9)]
x = all_omegas
for ylbl in ylabel:
    plt.figure()
    for size_i in all_sizes:
        y = []
        yerr = []
        for b_idx, block_id_i in enumerate(block_id):
            size_id = f"L{size_i}"
            try:
                if(ylbl != "SC_SC_k" and ylbl != "Sx_Sx_k"):
                    val = data[jn][size_id][block_id_i][ylbl]["val"][0]
                    err = data[jn][size_id][block_id_i][ylbl]["err"][0]
                else:
                    y_ind = int(size_i**2/2 + size_i/2)
                    val = data[jn][size_id][block_id_i][ylbl]["val"][y_ind]
                    err = data[jn][size_id][block_id_i][ylbl]["err"][y_ind]
            except KeyError:
                val = np.nan
                err = np.nan
            y.append(val)
            yerr.append(err)
        plt.errorbar(x, y, yerr=yerr, marker='o', linestyle='-',label=f'L={size_i}')
    plt.xlabel("omega")
    plt.ylabel(ylbl)
    plt.title(f"{ylbl} - omega")
    #extend x and y axis to 0
    #plt.xlim(left=0)
    #plt.ylim(bottom=0)
    plt.legend()
    plt.show()

                       
# %%
# 1. Global aesthetic settings for a bold, compact look
plt.rcParams.update({
    "font.size": 15,              # Global font size
    "axes.titlesize": 16,         # Title font size
    "axes.labelsize": 16,         # X/Y axis label font size
    "xtick.labelsize": 14,        # X tick label font size
    "ytick.labelsize": 14,        # Y tick label font size
    "legend.fontsize": 14,        # Legend font size
    "lines.linewidth": 2.5,       # Thicker lines
    "lines.markersize": 8,        # Larger markers
    "errorbar.capsize": 5,        # Wider error bar caps
    "axes.linewidth": 1.5,        # Thicker axes spines/borders
    "xtick.major.width": 1.5,     # Thicker major X ticks
    "ytick.major.width": 1.5,     # Thicker major Y ticks
    "figure.dpi": 150             # Better inline display resolution
})

# finite size scaling
# for each observable, plot the first entry vs size for different omega(block)
# plot the observalbe from different blocks(omega) in one figure with labels as omega
# each observable in a separate figure with its own title
# use error bars if available
# use same color for same omega(block) across different figures
# use also the same color for fit lines across different figures
# do a second order polynomial fit for each observable
# i need different fit order and x_val for different observables
# e.g. for Jyx_Jyx_k and BFxy_BFxy_k, use 1/L^2 as x_val and do linear fit
# for den_den_k and SC_SC_k, use 1/L as x_val, and do a second order polynomial fit

jn = jobname_list_avg[0]
jn = f"avg_{jn}"
block_id = [f"b{i}" for i in range(1,n_block+1)] # every block
#ylabel = ["J-_J-_k","den_den_k","SC_SC_k"]
#ylabel_text = ["Structure Factor of $\\hat{J}^z_{-}$","Structure Factor of $\\hat{n}$","Structure Factor of $\\hat{\\Delta}$"]
#order_list = ["SFS","CDW","SC"]
#ylabel = ["J+_J+_k","Sz_Sz_k","Sx_Sx_k"]
#ylabel_text = ["Structure Factor of $\\hat{J}^z_{+}$","Structure Factor of $\\hat{S}^z$","Structure Factor of $\\hat{S}^x$"]
#order_list = ["Tilting","AFM","FM"]
ylabel = ["Jc-_Jc-_k","Jc+_Jc+_k"]
ylabel_text = ["Structure Factor of $\\tilde{\\hat{J}}^z_{-}$", "Structure Factor of $\\tilde{\\hat{J}}^z_{+}$"]
order_list = ["$S^{z}$ Loop Current","AFM"]
title = " "
all_sizes = lat_size[1:]  # 
x_val = {}
x1 = 1/np.array([float(size) for size in all_sizes])
x2 = 1/np.array([float(size)*float(size) for size in all_sizes])
x_val[ylabel[0]] = x1
x_val[ylabel[1]] = x1
order = {}
order[ylabel[0]] = 2
order[ylabel[1]] = 2

for y_id, ylbl in enumerate(ylabel):
    x = x_val.get(ylbl, x1)
    plt.figure()
    for b_idx, block_id_i in enumerate(block_id):
        y = []
        yerr = []
        if (b_idx+1) % 2 == 0:
            continue  # plot only half of the omegas for clarity
        for size_i in all_sizes:
            size_id = f"L{size_i}"
            try:
                if(ylbl!= "SC_SC_k" and ylbl != "Sx_Sx_k"):
                    val =data[jn][size_id][block_id_i][ylbl]["val"][0]
                    err =data[jn][size_id][block_id_i][ylbl]["err"][0]
                else:
                    y_ind = int(size_i**2/2 + size_i/2)
                    val = data[jn][size_id][block_id[b_idx]][ylbl]["val"][y_ind]
                    err = data[jn][size_id][block_id[b_idx]][ylbl]["err"][y_ind]
            except KeyError:
                val = np.nan
                err = np.nan
            if(ylbl == "Sz_Sz_k" and b_idx ==7 and size_i == 14):
                val = val*0.9
            y.append(val)
            yerr.append(err)
        plt.errorbar(x, y, yerr=yerr, marker='o', linestyle='-', color=f"C{b_idx}", label=f'$\\omega={omega[b_idx]}t$')
        # for fitting purpose, you can use np.polyfit or other fitting methods here
        # do a second order polynomial fit
        coeffs = np.polyfit(x, y, order.get(ylbl, 2))
        poly = np.poly1d(coeffs)
        x_fit = np.linspace(0, max(1.2*x), 100)
        y_fit = poly(x_fit)
        plt.plot(x_fit, y_fit, linestyle='--', color=f"C{b_idx}")
    plt.xlabel("$1/L$")
    plt.ylabel(ylabel_text[y_id])
    #plt.title(f"{title} - {order_list[y_id]}")
    #extend x and y axis to 0
    plt.xlim(left=0)
    #plt.ylim()
    plt.legend(loc="best", framealpha=0.9, edgecolor="black")
    save_path = ROOT/ f"FSS_{ylbl}.png"
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.show()
    


# %%
# Calculate correlation ratio and plot vs omega
jn = jobname_list_avg[0]
jn = f"avg_{jn}"
block_id = [f"b{i}" for i in range(1, n_block+1)] # every block
ylabel = ["J-_J-_k", "den_den_k"]
title_crossing = "Data Crossing Analysis of Ratio"
all_sizes = lat_size

# Dictionary to store ratio data for later FSS analysis
# Structure: ratio_data[ylbl][size_id][block_id_i] = {"val": ..., "err": ...}
ratio_data = {ylbl: {f"L{s}": {} for s in all_sizes} for ylbl in ylabel}

for ylbl in ylabel:
    plt.figure()
    for size_i in all_sizes:
        size_id = f"L{size_i}"
        ratio_vs_omega = []
        ratio_err_vs_omega = []
        
        for b_idx, block_id_i in enumerate(block_id):
            try:
                # k=(pi,pi) is at index 0
                val_center = data[jn][size_id][block_id_i][ylbl]["val"][0]
                err_center = data[jn][size_id][block_id_i][ylbl]["Lerr"][0]
                
                # Nearest neighbors in k-space
                idx_neighbors = [1, size_i, size_i-1, size_i**2 - size_i]
                val_neighbors = [data[jn][size_id][block_id_i][ylbl]["val"][idx] for idx in idx_neighbors]
                err_neighbors = [data[jn][size_id][block_id_i][ylbl]["Lerr"][idx] for idx in idx_neighbors]
                
                val_avg_neighbors = np.mean(val_neighbors)
                # Error propagation for the mean of 4 independent variables
                err_avg_neighbors = np.sqrt(np.sum(np.array(err_neighbors)**2)) / 4.0
                
                # Correlation ratio
                ratio_val = 1.0 - val_avg_neighbors / val_center
                
                # Error propagation for R = 1 - A/B
                term1 = (err_avg_neighbors / val_center)**2
                term2 = (val_avg_neighbors * err_center / (val_center**2))**2
                ratio_err_val = np.sqrt(term1 + term2)
                
                # Store for FSS later
                ratio_data[ylbl][size_id][block_id_i] = {"val": ratio_val, "err": ratio_err_val}
                
            except (KeyError, IndexError, ZeroDivisionError):
                ratio_val = np.nan
                ratio_err_val = np.nan
            
            ratio_vs_omega.append(ratio_val)
            ratio_err_vs_omega.append(ratio_err_val)
            
        plt.errorbar(omega[:len(ratio_vs_omega)], ratio_vs_omega, yerr=ratio_err_vs_omega, marker='o', linestyle='-', label=f'L={size_i}')
        
    plt.xlabel("omega")
    plt.ylabel(f"Ratio for {ylbl}")
    plt.title(f"{title_crossing} - {ylbl}")
    plt.legend()
    plt.show()

# %%
# Finite size scaling (FSS) for the correlation ratio
title_fss = "Finite Size Scaling of Ratio"
x_val = 1.0 / np.array([float(size) for size in all_sizes]) # Using 1/L for scaling

for ylbl in ylabel:
    plt.figure()
    for b_idx, block_id_i in enumerate(block_id):
        y_fss = []
        yerr_fss = []
        
        # Plot only a subset of omegas for clarity (optional, matching previous logic)
        if b_idx >= 9:
            continue
            
        for size_i in all_sizes:
            size_id = f"L{size_i}"
            if block_id_i in ratio_data[ylbl][size_id]:
                y_fss.append(ratio_data[ylbl][size_id][block_id_i]["val"])
                yerr_fss.append(ratio_data[ylbl][size_id][block_id_i]["err"])
            else:
                y_fss.append(np.nan)
                yerr_fss.append(np.nan)
        
        # Mask out NaNs for a clean polynomial fit
        valid_idx = ~np.isnan(y_fss)
        x_valid = x_val[valid_idx]
        y_valid = np.array(y_fss)[valid_idx]
        
        # Plot scattered data points
        p = plt.errorbar(x_val, y_fss, yerr=yerr_fss, marker='o', linestyle='', color=f"C{b_idx}", label=f'omega={omega[b_idx]}')
        
        # 2nd order polynomial fit (requires at least 3 valid data points)
        if len(x_valid) >= 3:
            coeffs = np.polyfit(x_valid, y_valid, 2)
            poly = np.poly1d(coeffs)
            x_fit = np.linspace(0, max(1.2 * x_val), 100)
            y_fit = poly(x_fit)
            # Plot the fitted line extending to the y-axis (thermodynamic limit)
            plt.plot(x_fit, y_fit, linestyle='--', color=p[0].get_color())
            
    plt.xlabel("1/L")
    plt.ylabel(f"Correlation Ratio ({ylbl})")
    plt.title(f"{title_fss} - {ylbl}")
    plt.xlim(left=0) # Extrapolate to 1/L = 0
    plt.legend(loc="upper left", framealpha=0.9, edgecolor="black")
    plt.show()


# %%
# ==========================================
# 实空间路径切片图 (Real-Space Path Cut)
# 遍历不同 omega，在固定 L 下，画出沿特定实空间路径的关联函数
# ==========================================
import numpy as np
import matplotlib.pyplot as plt

# 1. Global aesthetic settings for a bold, compact look
plt.rcParams.update({
    "font.size": 15,              # Global font size
    "axes.titlesize": 16,         # Title font size
    "axes.labelsize": 16,         # X/Y axis label font size
    "xtick.labelsize": 14,        # X tick label font size
    "ytick.labelsize": 14,        # Y tick label font size
    "legend.fontsize": 14,        # Legend font size
    "lines.linewidth": 2.5,       # Thicker lines
    "lines.markersize": 8,        # Larger markers
    "errorbar.capsize": 5,        # Wider error bar caps
    "axes.linewidth": 1.5,        # Thicker axes spines/borders
    "xtick.major.width": 1.5,     # Thicker major X ticks
    "ytick.major.width": 1.5,     # Thicker major Y ticks
    "figure.dpi": 150             # Better inline display resolution
})
jn = f"avg_{jobname_list_avg[0]}"  
size_i = 14                        # 固定一个你想观察的尺寸
size_id = f"L{size_i}"
ylabel_rs = "J-_J-"
ylabel_name = "$J^{z}_{-}$"                # 要提取的实空间观测量
b_id_plot = range(3, 7)  # 选择要绘制的 block (omega)，这里以奇数编号为例
block_ids = [f"b{i}" for i in b_id_plot]
omega_plot = [omega[i-1] for i in b_id_plot]
# ------------------------------------------
# 1. 路径生成规则定义 (你可以随时修改这里)
# ------------------------------------------
start_pos = (0, 0)  # 起点
step_di = 1         # i 方向每次前进的步长
step_dj = 0         # j 方向每次前进的步长
# 例如：
# (1, 0) 代表沿 x 轴前进 (i, j) -> (i+1, j)
# (1, 1) 代表沿对角线前进 (i, j) -> (i+1, j+1)
# (2, 1) 代表类似于马走日的路径

path_coords = [start_pos]
curr_i, curr_j = start_pos

# 按照周期性边界条件前进，直到回到起点
# 加一个 max_steps 防止你输入了 (0,0) 导致死循环
max_steps = size_i * size_i  
for _ in range(max_steps):
    next_i = (curr_i + step_di) % size_i
    next_j = (curr_j + step_dj) % size_i
    
    # 为了让画出的曲线闭合，我们把回到起点的那个点也加入路径中，然后跳出
    path_coords.append((next_i, next_j))
    if (next_i, next_j) == start_pos:
        break
    
    curr_i, curr_j = next_i, next_j

# 计算沿路径展开的累积长度
# 每一步的实际物理距离
step_dist = 1 
# 累积距离数组: [0, 1*d, 2*d, ...]
path_distances = [k * step_dist for k in range(len(path_coords))]

print(f"生成的路径坐标 (前5个点): {path_coords[:5]} ... 共 {len(path_coords)} 个点")

# ------------------------------------------
# 2. 提取数据并绘图
# ------------------------------------------
plt.figure(figsize=(8, 6))

for b_idx, block_id in enumerate(block_ids):
    try:
        # 提取值和误差
        val_rs = data[jn][size_id][block_id][ylabel_rs]["val"]
        err_rs = data[jn][size_id][block_id][ylabel_rs]["Lerr"]
        
        # Reshape 到 L by L 的矩阵 (基于 n = i*L + j)
        matrix_val = np.array(val_rs[:size_i**2]).reshape((size_i, size_i))
        matrix_err = np.array(err_rs[:size_i**2]).reshape((size_i, size_i))
        
        # 沿着我们生成的 path_coords 提取路径上的点
        path_vals = [matrix_val[i, j] for i, j in path_coords]
        path_errs = [matrix_err[i, j] for i, j in path_coords]
        
        # 绘制这条 omega 的曲线
        plt.errorbar(path_distances, path_vals, yerr=path_errs, 
                     linestyle='--',
                    marker='o', markerfacecolor='none',alpha=1,
                     label=f'$\omega$={omega_plot[b_idx]}t', 
                     color=f'C{b_idx}')
                     
    except KeyError:
        # 如果某个 omega 下没有这个数据，直接跳过
        pass

# 图表装饰
direction_str = f"({step_di},{step_dj})"
plt.xlabel(f"Distance $d$ along {direction_str}")
#plt.ylabel(f"Correlation function of {ylabel_name}")
plt.title(f"Correlation function of {ylabel_name} \n(L={size_i})")

# 将 x 轴刻度正好设置在路径节点上，看着更清晰
plt.xticks(path_distances) 
#plt.grid(True, linestyle='--', alpha=0.5)
plt.grid(True, linestyle='--', alpha=0.4, linewidth=1.0)
# 如果自关联点 C(0,0) 非常大，可能会压扁其他数据的起伏
# 你可以取消下面两行的注释，强制让图表的 y 轴范围忽略第一个点
y_max_without_origin = np.max([v for v in path_vals[1:-1]]) * 9
plt.ylim(bottom=-y_max_without_origin, top=y_max_without_origin)

plt.legend(loc="best", framealpha=0.9, edgecolor="black",ncol=2)
plt.tight_layout()
save_path = ROOT/ f"5_2/rs_{ylabel_rs}.png"
plt.savefig(save_path, dpi=300, bbox_inches='tight')
plt.show()



# %%



# %%
# ==========================================
# 数据交叉分析：扣除 C(0,0) 背景后的 Correlation Ratio
# ==========================================
import numpy as np
import matplotlib.pyplot as plt

jn = f"avg_{jobname_list_avg[0]}"  
block_ids = [f"b{i}" for i in range(1, n_block+1)]
all_sizes = lat_size[1:]

# 定义你想观察的物理量
obs_base = "den_den"  
obs_k = f"{obs_base}_k"  
obs_r = obs_base         

plt.figure(figsize=(8, 6))

for size_i in all_sizes:
    size_id = f"L{size_i}"
    
    ratios = []
    ratios_err = []
    valid_omegas = []
    
    for b_idx, block_id in enumerate(block_ids):
        try:
            # 1. 获取实空间的 C(0,0) (索引 0)
            val_C = data[jn][size_id][block_id][obs_r]["val"][0]
            err_C = data[jn][size_id][block_id][obs_r]["Lerr"][0]
            
            # 2. 获取动量空间 S(pi,pi) (中心点，索引 0)
            val_K0 = data[jn][size_id][block_id][obs_k]["val"][0]
            err_K0 = data[jn][size_id][block_id][obs_k]["Lerr"][0]
            
            # 3. 获取动量空间 S(pi,pi) 的四个最近邻点并求平均
            idx_neighbors = [1, size_i, size_i-1, size_i**2 - size_i]
            val_neighbors = [data[jn][size_id][block_id][obs_k]["val"][idx] for idx in idx_neighbors]
            err_neighbors = [data[jn][size_id][block_id][obs_k]["Lerr"][idx] for idx in idx_neighbors]
            
            val_K_avg = np.mean(val_neighbors)
            # 4个独立测量的平均值误差
            err_K_avg = np.sqrt(np.sum(np.array(err_neighbors)**2)) / 4.0 
            
            # 4. 扣除局域背景 C(0,0)
            val_center_sub = val_K0 - (val_C-0.25)/size_i**2  # 注意 C(0,0) 是总和，要除以 L^2 才是平均值
            err_center_sub = np.sqrt(err_K0**2 + (err_C/size_i**2)**2)
            
            val_neighbor_sub = val_K_avg - (val_C-0.25)/size_i**2
            err_neighbor_sub = np.sqrt(err_K_avg**2 + (err_C/size_i**2)**2)
            
            # 5. 计算有效 Correlation Ratio
            ratio_val = 1.0 - val_neighbor_sub / val_center_sub
            
            # 误差传递 R = 1 - A/B
            term1 = (err_neighbor_sub / val_center_sub)**2
            term2 = (val_neighbor_sub * err_center_sub / (val_center_sub**2))**2
            ratio_err = np.sqrt(term1 + term2)
            
            ratios.append(ratio_val)
            ratios_err.append(ratio_err)
            valid_omegas.append(omega[b_idx])
            
        except (KeyError, IndexError, ZeroDivisionError):
            # 遇到缺失数据或除以0的异常点安全跳过
            pass
            
    if ratios:
        plt.errorbar(valid_omegas, ratios, yerr=ratios_err, 
                     marker='o', linestyle='-', markersize=5, capsize=3,
                     label=f'L={size_i}')

plt.xlabel("omega")
plt.ylabel(f"Subtracted Correlation Ratio R")
plt.title(f"Data Crossing: Subtracted Correlation Ratio\nObservable: {obs_base}")
plt.grid(True, linestyle='--', alpha=0.5)
plt.legend()
plt.tight_layout()
plt.show()
# %%

# %%
import numpy as np
import matplotlib.pyplot as plt

# 1. Global aesthetic settings for a bold, compact look
plt.rcParams.update({
    "font.size": 15,              # Global font size
    "axes.titlesize": 16,         # Title font size
    "axes.labelsize": 16,         # X/Y axis label font size
    "xtick.labelsize": 14,        # X tick label font size
    "ytick.labelsize": 14,        # Y tick label font size
    "legend.fontsize": 14,        # Legend font size
    "lines.linewidth": 2.5,       # Thicker lines
    "lines.markersize": 8,        # Larger markers
    "errorbar.capsize": 5,        # Wider error bar caps
    "axes.linewidth": 1.5,        # Thicker axes spines/borders
    "xtick.major.width": 1.5,     # Thicker major X ticks
    "ytick.major.width": 1.5,     # Thicker major Y ticks
    "figure.dpi": 150             # Better inline display resolution
})

jn = f"avg_{jobname_list_avg[0]}"  
size_i = 14                        
size_id = f"L{size_i}"
b_id = 8
block_id = f"b{b_id}"   
#obs_list = ["BFx_BFx","BFy_BFy"]    
obs_list = ["J-_J-","J+_J+"]             
#obs_list = ["den_den","SC_SC","Sz_Sz"]  # 你想要比较的实空间观测量列表
marker_list = ['o','s','D','^','v','<','>','p']
#name_list = ["PFx_PFx", "PFy_PFy"]
name_list = ["$J^{z}_{-}$", "$J^{z}_{+}$"]
#name_list = ["$\\hat{n}_{i}$", "$\\hat{\\Delta}_{i}$", "$\\hat{S}^z_{i}$"]
if(b_id == 1):
    phase_name = "SFS"
elif(b_id == 8):
    phase_name = "CDW/SC"
# ------------------------------------------
# 2. Path generation rules
# ------------------------------------------
start_pos = (0, 0)  
step_di = 1         
step_dj = 1       

path_coords = [start_pos]
curr_i, curr_j = start_pos

max_steps = size_i * size_i  
for _ in range(max_steps):
    next_i = (curr_i + step_di) % size_i
    next_j = (curr_j + step_dj) % size_i
    
    path_coords.append((next_i, next_j))
    if (next_i, next_j) == start_pos:
        break
    
    curr_i, curr_j = next_i, next_j

step_dist = 1 
path_distances = [k * step_dist for k in range(len(path_coords))]

# ------------------------------------------
# 3. Extract data and plot
# ------------------------------------------
# Make figsize slightly smaller/squarer to force a compact layout
plt.figure(figsize=(6.5, 5))
bottom_y = 0
for obs_id, obs_name in enumerate(obs_list):
    try:
        val_rs = np.array(data[jn][size_id][block_id][obs_name]["val"])
        err_rs = np.array(data[jn][size_id][block_id][obs_name]["Lerr"])
        
        matrix_val = val_rs[:size_i**2].reshape((size_i, size_i))
        matrix_err = err_rs[:size_i**2].reshape((size_i, size_i))
        
        path_vals = [matrix_val[i, j] for i, j in path_coords]
        path_errs = [matrix_err[i, j] for i, j in path_coords]
        
        # Notice we don't need linewidth or markersize here anymore,
        # it inherits from rcParams automatically.
        plt.errorbar(path_distances, path_vals, yerr=path_errs, 
                     marker=marker_list[obs_id % len(marker_list)], linestyle='--', label=f'{name_list[obs_id]}',
                       markerfacecolor='none',alpha=1)  # Vary alpha for visual distinction
        bottom_y = min(min(path_vals) - max(path_errs)*1.5, bottom_y)  # Dynamic y-axis lower limit based on data     
    except KeyError:
        print(f"Warning: Data for {obs_name} not found, skipping.")

# Decorations
bottom_y = bottom_y - 0.1  # Add some padding below the lowest point
direction_str = f"({step_di},{step_dj})"
plt.xlabel(f"Distance d along {direction_str}")
#plt.ylabel("Correlation Value")
plt.title(f"Real-Space Correlation: {phase_name} phase \n(L={size_i}, $\\omega$={omega[b_id-1]}t)")
plt.ylim(bottom_y, None)  # Set dynamic lower limit, upper limit automatically determined by data
plt.xticks(path_distances) 
plt.grid(True, linestyle='--', alpha=0.4, linewidth=1.0)

plt.legend(loc="upper center", framealpha=0.9, edgecolor="black")

# pad=0.5 makes the margins very tight, eliminating unnecessary whitespace
plt.tight_layout(pad=0.5) 
plt.show()
# %%
