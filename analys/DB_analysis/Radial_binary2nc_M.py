#%%
import numpy as np
import pandas as pd
import glob
import os 
import xarray as xr
import sys

# %%
idx = int (sys.argv[1])
heig = str(np.arange(5, 210, 10)[idx]).zfill(3)
file_dir = os.path.dirname(os.path.realpath(__file__))
files_U = sorted(glob.glob(file_dir[:-6] + "u_slice/U%st=*"%heig))
files_V = sorted(glob.glob(file_dir[:-6] + "u_slice/V%st=*"%heig))
files_W = sorted(glob.glob(file_dir[:-6]
 + "u_slice/W%st=*"%heig))
vaR = "M"
# %% save the variable
dir_out = file_dir + "/Radial/"
if not os.path.exists(dir_out):
    os.makedirs(dir_out)
# %% the array now is t, y, z
xl = np.arange(2049)
zl = np.arange(2049)
tt = pd.to_timedelta(np.arange(0, 3000, 4), unit="s")

ll_bins = np.arange(0, 500, 2)
ll_bins_label = np.arange(1, 499, 2)
#%%
for i in range(0, len(files_U)):
    print(i)
    sqe = str(i).zfill(2)
    USlice = np.fromfile(files_U[i], dtype = float).reshape(25, 2049, 2049)
    VSlice = np.fromfile(files_V[i], dtype = float).reshape(25, 2049, 2049)
    WSlice = np.fromfile(files_W[i], dtype = float).reshape(25, 2049, 2049)
    MSlice = (USlice**2 + VSlice**2 + WSlice**2)**0.5
    b_xr = xr.DataArray(data=MSlice, dims=["time", "x", "z"], coords=[tt[i*25:25*(i+1)], xl, zl])
    b_xr = b_xr.assign_coords(ll = ((b_xr.x - 500)**2 + (b_xr.z - 1024)**2)**0.5)
    # b_xr = b_xr.assign_coords(ang = np.arctan2(b_xr.x - 500, b_xr.z - 1024)* 180 / np.pi)
    b_xr_RA = b_xr.groupby_bins("ll", ll_bins, labels=ll_bins_label).mean()
    b_xr_RA = b_xr_RA.dropna("ll_bins")
    b_ds_RA = b_xr_RA.to_dataset(name = "tmp")
    # np.save(dir_out + "TX_%s_%s.npy"%(vaR, sqe), Slice)
    b_ds_RA.to_netcdf(dir_out + "%sRA_y%s_t%s.nc"%(vaR, heig, sqe))