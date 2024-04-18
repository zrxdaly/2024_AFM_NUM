#%% 
import numpy as np
import pandas as pd
import xarray as xr
import glob
import os
import sys
#%%
file_dir = os.path.dirname(os.path.realpath(__file__))
b_AG_file = sorted(glob.glob(file_dir + "/Radial_profile/Radial_Bprof.nc"))
b_AG_y = xr.load_dataset(b_AG_file[0], engine = "netcdf4")
#%%
rb_AG_REF = b_AG_y.tmp[:,:150,:]
TA_rb_AG_REF = rb_AG_REF.mean("time")
#%% save the variables for plotting
np.save(file_dir + "/Radial_profile/Rb_AG_RE.npy", TA_rb_AG_REF)
#%%
# PAV_rb_AG_OP = b_AG_y.tmp.values[:,150+72*0:150+72*(0+1), :]
# for j in range(1, 4):
#     print(j)
#     PAV_rb_AG_OP = PAV_rb_AG_OP + b_AG_y.tmp.values[:,150+72*j:150+72*(j+1), :]
# PAV_rb_AG_OP = PAV_rb_AG_OP/4
PAV_rb_AG_OP = b_AG_y.tmp.values[:,150:, :]
# %% ---------------------------------------------------------- #
# -- AG* zone averaging at whole period of time # 
# % ----------------------------------------------------------- #
Rb_AG_OP = PAV_rb_AG_OP.mean(axis = 1)

#%% save the variables for plotting
np.save(file_dir + "/Radial_profile/Rb_AG_OP.npy", Rb_AG_OP)

# %%
