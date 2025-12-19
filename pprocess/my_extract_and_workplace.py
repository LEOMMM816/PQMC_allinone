

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
jobname_list = ["O1623","BCO1623"]
n_block = 8
# omega in 1.6,1.7...2.3, kept to 2 significant digits
omega = [round(1.6 + 0.1*i,2) for i in range(n_block)]
lat_size = [4, 6, 8, 10,12]
data = {}
for jn in jobname_list:
    data[jn] = data.get(jn, {})
    for size in lat_size:
        file_name = f"{jn}L{size}"
        size_name = f"L{size}"
        data[jn][size_name] = data[jn].get(size_name, {})
        base_path = ROOT/f"./{file_name}/data/out_files/"
        if not base_path.exists():
            print(f"{base_path} doesn't exist, skipped")
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
                data[jn][size_name][block_name][name] = data[jn][size_name][block_name].get(name, {"vals": [], "errs": []})
                try:
                    val = float(parts[2])
                    err = float(parts[3])
                    data[jn][size_name][block_name][name]["vals"].append(val)
                    data[jn][size_name][block_name][name]["errs"].append(err)
                except ValueError:
                    if name != "Name":
                        print(f"Could not find value for {name} in {jn} {size_name} {block_name}, skipping...")
                    continue
   

# %% 
# average the data from with BC and without BC
jobname_list_avg = ["O1623"]
data_avg = {}
for jn in jobname_list_avg:
    jn1 = jn  # without BC
    jn2 = "BC" + jn  # with BC
    data_avg[jn] = {}
    for size in lat_size:
        size_name = f"L{size}"
        data_avg[jn][size_name] = {}
        for block in range(1, n_block+1):
            block_name = f"b{block}"
            data_avg[jn][size_name][block_name] = {}
            keys = set(data[jn1][size_name][block_name].keys()).intersection(set(data[jn2][size_name][block_name].keys()))
            for key in keys:
                vals1 = np.array(data[jn1][size_name][block_name][key]["vals"])
                errs1 = np.array(data[jn1][size_name][block_name][key]["errs"])
                vals2 = np.array(data[jn2][size_name][block_name][key]["vals"])
                errs2 = np.array(data[jn2][size_name][block_name][key]["errs"])
                avg_vals = (vals1 + vals2) / 2
                avg_errs = np.sqrt(errs1**2 + errs2**2) / 2
                data_avg[jn][size_name][block_name][key] = {"vals": avg_vals.tolist(), "errs": avg_errs.tolist()}
#plot some data from data_avg for verification
example_name = "Jyx_Jyx_k"
jn1 = jn  # without BC
jn2 = "BC" + jn  # with BC
for size in lat_size:
    size_name = f"L{size}"
    block_name = "b4"
    vals1 = np.array(data[jn1][size_name][block_name][example_name]["vals"])
    errs1 = np.array(data[jn1][size_name][block_name][example_name]["errs"])
    vals2 = np.array(data[jn2][size_name][block_name][example_name]["vals"])
    errs2 = np.array(data[jn2][size_name][block_name][example_name]["errs"])
    plt.errorbar(size,vals1[0], yerr=errs1[0], marker='o', linestyle='-', color='blue', label=f'{jn1} {size_name} {block_name} {example_name}')
    plt.errorbar(size,vals2[0], yerr=errs2[0], marker='o', linestyle='-', color='red', label=f'{jn2} {size_name} {block_name} {example_name}')
    
# %%
# plot certain observable vs real(momentum) space for a given lattice size and block(omega)
jn = jobname_list_avg[0]
size_i = 12
size_id = f"L{size_i}"
block_id = "b1"
xlabel = "x"
#ylabel = "BFxy_BFxy_k"
# ylabel = "SC_SC_k"
ylabel = "Jxy_Jxy_k"
title = "Test plot"
x = range(size_i**2)  # example x values
y = data[jn][size_id][block_id][ylabel]["vals"]
y[int(size_i**2/2 + size_i/2)] = 0  # set the middle point to 0 for better visualization
yerr = data[jn][size_id][block_id][ylabel]["errs"]
plt.errorbar(x, y, yerr=yerr, marker='o', linestyle='-', label=f'{size_id} {block_id} {ylabel}')
plt.xlabel(xlabel)
plt.ylabel(ylabel)
plt.title(title)
plt.legend()
plt.show()


# %%
# plot observable vs different blocks(omega) for lattice sizes
jn = jobname_list[1]
block_id = [f"b{i}" for i in range(1,9)]
ylabel = ["Jxy_Jxy_k", "den_den_k", "SC_SC_k"]
title = "Finite Size Scaling Analysis"
all_sizes = lat_size
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
                if(ylbl != "SC_SC_k"):
                    val = data[jn][size_id][block_id_i][ylbl]["vals"][0]
                    err = data[jn][size_id][block_id_i][ylbl]["errs"][0]
                else:
                    y_ind = int(size_i**2/2 + size_i/2)
                    val = data[jn][size_id][block_id_i][ylbl]["vals"][y_ind]
                    err = data[jn][size_id][block_id_i][ylbl]["errs"][y_ind]
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
    plt.ylim(bottom=0)
    plt.legend()
    plt.show()

                       
# %%
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
block_id = [f"b{i}" for i in range(1,n_block+1)] # every block
ylabel = ["Jxy_Jxy_k","den_den_k"]
title = "Finite Size Scaling Analysis"
all_sizes = lat_size  # 
x_val = {}
x1 = 1/np.array([float(size) for size in all_sizes])
x2 = 1/np.array([float(size)*float(size) for size in all_sizes])
x_val[ylabel[0]] = x2
x_val[ylabel[1]] = x1
order = {}
order[ylabel[0]] = 1
order[ylabel[1]] = 2

for ylbl in ylabel:
    x = x_val[ylbl]
    plt.figure()
    for b_idx, block_id_i in enumerate(block_id):
        y = []
        yerr = []
        if b_idx % 2 ==0:
            continue  # plot only half of the omegas for clarity
        for size_i in all_sizes:
            size_id = f"L{size_i}"
            try:
                if(ylbl!= "SC_SC_k"):
                    val = data[jn][size_id][block_id_i][ylbl]["vals"][0]
                    err = data[jn][size_id][block_id_i][ylbl]["errs"][0]
                else:
                    y_ind = int(size_i**2/2 + size_i/2)
                    val = data[jn][size_id][block_id_i][ylbl]["vals"][y_ind]
                    err = data[jn][size_id][block_id_i][ylbl]["errs"][y_ind]
            except KeyError:
                val = np.nan
                err = np.nan
            y.append(val)
            yerr.append(err)
        plt.errorbar(x, y, yerr=yerr, marker='o', linestyle='-', color=f"C{b_idx}", label=f'omega={omega[b_idx]}')
        # for fitting purpose, you can use np.polyfit or other fitting methods here
        # do a second order polynomial fit
        coeffs = np.polyfit(x, y, order[ylbl])
        poly = np.poly1d(coeffs)
        x_fit = np.linspace(0, max(1.2*x), 100)
        y_fit = poly(x_fit)
        plt.plot(x_fit, y_fit, linestyle='--', color=f"C{b_idx}")
    plt.xlabel("Lattice Size")
    plt.ylabel(ylbl)
    plt.title(f"{title} - {ylbl}")
    #extend x and y axis to 0
    plt.xlim(left=0)
    #plt.ylim()
    plt.legend()
    plt.show()
    
    
  # %%

