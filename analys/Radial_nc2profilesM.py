#%%
import numpy as np
import pandas as pd
import glob
import os 
import xarray as xr
import sys
# %% save the variable
file_dir = os.path.dirname(os.path.realpath(__file__))
dir_out = file_dir + "/Radial_profile/"
if not os.path.exists(dir_out):
    os.makedirs(dir_out)
# %%
TT = str(np.arange(0, 30, 1)[0]).zfill(2)
files_nc = sorted(glob.glob(file_dir + "/Radial/MRA_y*_t%s.nc"%TT))
B_prof = xr.load_dataset(files_nc[0])
for i in range(1, len(files_nc)):
    B0 = xr.load_dataset(files_nc[i])
    B_prof = xr.concat([B_prof, B0], 'y')
yy = np.arange(0.5, 21., 1.)
B_prof = B_prof.assign_coords(y = yy)

for idx in range(1, 30):
    TT = str(np.arange(0, 30, 1)[idx]).zfill(2)
    files_nc = sorted(glob.glob(file_dir + "/Radial/MRA_y*_t%s.nc"%TT))
    B_prof0 = xr.load_dataset(files_nc[0])
    for i in range(1, len(files_nc)):
        B0 = xr.load_dataset(files_nc[i])
        B_prof0 = xr.concat([B_prof0, B0], 'y')
    B_prof0 = B_prof0.assign_coords(y = yy)
    ## after here we stack time together
    B_prof = xr.concat([B_prof, B_prof0], 'time')
#%%
B_prof.to_netcdf(dir_out + "Radial_Mprof.nc")
#%%












# %% the array now is t, y, z
xl = np.arange(2049)
zl = np.arange(2049)
tt = pd.to_timedelta(np.arange(3000), unit="s")

ll_bins = np.arange(0, 500, 2)
ll_bins_label = np.arange(1, 499, 2)
#%%
for i in range(0, len(files_Sl)):
    print(i)
    sqe = str(i).zfill(2)
    Slice = np.fromfile(files_Sl[i], dtype = float).reshape(25, 2049, 2049)/9.81*273.15
    b_xr = xr.DataArray(data=Slice, dims=["time", "x", "z"], coords=[tt[0*25:25*(0+1)], xl, zl])
    b_xr = b_xr.assign_coords(ll = ((b_xr.x - 500)**2 + (b_xr.z - 1024)**2)**0.5)
    # b_xr = b_xr.assign_coords(ang = np.arctan2(b_xr.x - 500, b_xr.z - 1024)* 180 / np.pi)
    b_xr_RA = b_xr.groupby_bins("ll", ll_bins, labels=ll_bins_label).mean()
    b_xr_RA = b_xr_RA.dropna("ll_bins")
    b_ds_RA = b_xr_RA.to_dataset(name = "tmp")
    np.save(dir_out + "TX_%s_%s.npy"%(vaR, sqe), Slice)
    b_ds_RA.to_netcdf(dir_out + "%sRA_y%s_t%s.nc"%(vaR, heig, sqe))