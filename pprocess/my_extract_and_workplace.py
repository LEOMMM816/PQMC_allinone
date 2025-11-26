

# %%
from pathlib import Path
import os

if "__file__" in globals():
    # 在正常脚本 / Run Python File / Run Current File in Interactive Window 时
    ROOT = Path(__file__).resolve().parent.parent
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
jobname = "epzO1"
n_block = 8
omega = [1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7]
lat_size = [4, 8, 10]
data = {}
for size in lat_size:
    file_name = f"{jobname}L{size}"
    run_name = f"L{size}"
    data[run_name] = {}
    base_path = ROOT / f"./results/{file_name}/data/out_files/"
    print(f"Loading data from {file_name} ...")
    for block in range(1, n_block+1):
        block_name = f"b{block}"
        data[run_name][block_name] = {}
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
            data[run_name][block_name][name] = data[run_name][block_name].get(name, {"vals": [], "errs": []})
            try:
                val = float(parts[2])
                err = float(parts[3])
                data[run_name][block_name][name]["vals"].append(val)
                data[run_name][block_name][name]["errs"].append(err)
            except ValueError:
                continue
# %%      


# %%
# Now data is loaded, we can plot
xlabel = "x"
ylabel = "Jxy_Jxy"
title = "Test plot"
size_id = "L10"
block_id = "b8"
# example plot for L4, block 1, observable Jxy_Jxy_k
x = range(100)  # example x values
y = data[size_id][block_id][ylabel]["vals"]
yerr = data[size_id][block_id][ylabel]["errs"]
plt.errorbar(x, y, yerr=yerr, marker='o', linestyle='-', label=f'{size_id} {block_id} {ylabel}')
plt.xlabel(xlabel)
plt.ylabel(ylabel)
plt.title(title)
plt.legend()
plt.show()
# %%