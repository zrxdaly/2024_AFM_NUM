#%% this script we compare the b2 with sim and exp
import numpy as np
from glob import glob
import os
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pandas as pd
import xarray as xr

plt.rc('font', family='serif')
plt.rc('xtick', labelsize='20')
plt.rc('ytick', labelsize='20')
plt.rc('text', usetex=True)

import scienceplots
plt.style.use(['science'])

# %% dir for different cases
base_dir = "/home/dai/Documents/Reserach/Fruit/Krabbendijke_experiment_2021_0507/VALIDATION_canopy/"
base_dir = os.path.dirname(os.path.realpath(__file__))
dirs_OP = sorted(glob(base_dir + "/W12/b2.npy"))
#%%
b2_OP = np.load(dirs_OP[0])
T2_OP = b2_OP[150:,:]
T2_REF = b2_OP[:150,:]
#%% prepare data from exp
in_folder = '/home/dai/Documents/Reserach/Fruit/Krabbendijke_experiment_2021_0507/DTS/calibration/'
b1_T_roll = xr.open_dataset(in_folder + "/2_contour_new/T_EF_save/" + "b1_T_roll.nc", engine = "netcdf4").rename({'__xarray_dataarray_variable__': 'tmpw'})
b2_T_roll = xr.open_dataset(in_folder + "/2_contour_new/T_EF_save/" + "b2_T_roll.nc", engine = "netcdf4").rename({'__xarray_dataarray_variable__': 'tmpw'})
# define later, could be different than exp 
OP_T_section = slice("2021-05-07T23:40:00", "2021-05-08T00:20:12")
RE_T_section = slice("2021-05-07T22:20:00", "2021-05-07T23:00:12")
b2_T_roll_OP  = b2_T_roll.sel(time =  OP_T_section)
b2_T_roll_RE  = b2_T_roll.sel(time =  RE_T_section)
b2_T_roll_OP_5S = b2_T_roll_OP.resample(time="4S").interpolate("linear").dropna(dim='time')
b2_T_roll_RE_5S = b2_T_roll_RE.resample(time="4S").interpolate("linear").dropna(dim='time')
# b1_EF.resample(time="1S").interpolate("linear").dropna(dim='time')

# %% time average b2 temperature 
y_exp = np.linspace(0,270,522)
x_exp = np.linspace(-150, 100, 30)
X_exp, Y_exp = np.meshgrid(x_exp, y_exp) 

y_sim = np.linspace(0, 270, 271)
x_sim = np.linspace(-150, 100, 251)
X_sim, Y_sim = np.meshgrid(x_sim, y_sim) 
#%%
fig, (ax1, ax2) = plt.subplots(1, 2, figsize = (12, 6), sharey=True)
# exp
DB_exp = np.mean(b2_T_roll_OP.tmpw, axis = 0) - np.mean(b2_T_roll_RE.tmpw,axis = 0)
con1 =ax1.contourf(X_exp, Y_exp, DB_exp, levels = np.arange(0, 6, 0.05), extend = "max", cmap = "gnuplot2")
tpp1 = ax1.contour(X_exp, Y_exp, DB_exp, levels = np.arange(0, 6, 2), colors = 'k', alpha = 0.90, linestyles = "solid")
ax1.clabel(tpp1, fmt='%2.1f', colors='w', fontsize=14)
# sim
DB_sim = np.mean(T2_OP, axis = 0) - np.mean(T2_REF, axis = 0)
DB_sim = DB_sim[5:276,5:256]
con2 =ax2.contourf(X_sim, Y_sim, DB_sim, levels = np.arange(0, 6, 0.05), extend = "max", cmap = "gnuplot2")
tpp1 = ax2.contour(X_sim, Y_sim, DB_sim, levels = np.arange(0, 6, 2), colors = 'k', alpha = 0.90, linestyles = "solid")
ax2.clabel(tpp1, fmt='%2.1f', colors='w', fontsize=14)
for c in con1.collections:
    c.set_edgecolor("face")
# for c in con2.collections:
#     c.set_edgecolor("face")
divider = make_axes_locatable(ax2)
cax = divider.append_axes("right", size="5%", pad=0.05)
hh = plt.colorbar(con1, cax=cax)
hh.ax.set_ylabel(r'$\Delta T$ [°C]',fontsize = 18)
plt.tight_layout()
# plt.savefig("b2_sim_num.png")
# %% create movie to diagnose if the wind direction is wrong or sth
# TimeF = b2_T_roll_OP_5S.index
def mp4_B2(DB_exp, DB_sim):
    dir_out = base_dir + "/b2_movies/"
    if not os.path.exists(dir_out):
            os.makedirs(dir_out)
    y_exp = np.linspace(0,270,522)
    x_exp = np.linspace(-150, 100, 30)
    X_exp, Y_exp = np.meshgrid(x_exp, y_exp) 

    y_sim = np.linspace(0, 270, 271)
    x_sim = np.linspace(-150, 100, 251)
    X_sim, Y_sim = np.meshgrid(x_sim, y_sim) 
    # exp
    for i in range(np.shape(DB_exp)[0]):
        with plt.ioff():
            timestring = str(b2_T_roll_OP_5S.time[i])[-19:].replace(':', '-')
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize = (12, 6), sharey=True)
            con1 =ax1.contourf(X_exp, Y_exp, DB_exp[i,:], levels = np.arange(0, 6, 0.05), extend = "max", cmap = "gnuplot2")
            tpp1 = ax1.contour(X_exp, Y_exp, DB_exp[i,:], levels = np.arange(0, 6, 2), colors = 'k', alpha = 0.90, linestyles = "solid")
            ax1.clabel(tpp1, fmt='%2.1f', colors='w', fontsize=14)
            #sim
            con2 =ax2.contourf(X_sim, Y_sim, DB_sim[i,:], levels = np.arange(0, 6, 0.05), extend = "max", cmap = "gnuplot2")
            tpp1 = ax2.contour(X_sim, Y_sim, DB_sim[i,:], levels = np.arange(0, 6, 2), colors = 'k', alpha = 0.90, linestyles = "solid")
            ax2.clabel(tpp1, fmt='%2.1f', colors='w', fontsize=14)
            for c in con1.collections:
                c.set_edgecolor("face")
            for c in con2.collections:
                c.set_edgecolor("face")
            divider = make_axes_locatable(ax2)
            cax = divider.append_axes("right", size="5%", pad=0.05)
            hh = plt.colorbar(con1, cax=cax)
            hh.ax.set_ylabel(r'$\Delta T$ [°C]',fontsize = 18)
            plt.tight_layout()
            plt.savefig(dir_out + "%s.png"%timestring)
            # plt.close()
        # break


# %%
DB_exp = b2_T_roll_OP_5S.tmpw.values - b2_T_roll_RE_5S.tmpw.values
T2_REF_AV = np.mean(T2_REF, axis = 0)
DB_sim = T2_OP- T2_REF_AV[np.newaxis, :,:]
DB_sim = DB_sim[:,5:276,5:256]
#%%
mp4_B2(DB_exp, DB_sim)
# %%
