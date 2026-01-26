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

# %%
basic_info = read_namelist('../data/data/reading_guide.nml', 'basic')
print(basic_info.keys())
obs_info= read_obs_info('../data/data/reading_guide.nml')
print(obs_info['BFx_BFy_k'])
# %%
with open('../data/data/out_core0.dat', 'r') as f:
    lines = f.readlines()
    # plot "BFx_BFy_k"[0] vs bin index
    obs_name = 'Jxy_Jxy_k'
    info = obs_info[obs_name]
    x = np.arange(basic_info['n_bins'])
    line_data = lines[info['offset_lo']-1]
    y = np.array([float(val) for val in line_data.split()])
    plt.plot(x, y)
    plt.xlabel('Bin index')
    plt.ylabel(obs_name + ' first k-point')
# %%
