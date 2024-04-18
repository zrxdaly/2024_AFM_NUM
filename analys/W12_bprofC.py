# %% ---------------------------------------------------------- #
# Here we compare the TH between sim and exp
# % ----------------------------------------------------------- #
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import datetime
from scipy.interpolate import InterpolatedUnivariateSpline
from matplotlib.dates import DateFormatter
import os

plt.rc('font', family='serif')
plt.rc('xtick', labelsize='20')
plt.rc('ytick', labelsize='20')
# plt.rc('text', usetex=True)

# %%
in_folder = '/home/dai/Documents/Reserach/Fruit/Krabbendijke_experiment_2021_0507/DTS/calibration'
#%% get the OP profiles
E06_OP_xr = xr.open_dataset(in_folder + "/3_Tower_new/Tower_xr/" + "E06_OP_xr.nc", engine = "netcdf4").rename({'__xarray_dataarray_variable__': 'tmpw'})
W12_OP_xr = xr.open_dataset(in_folder + "/3_Tower_new/Tower_xr/" + "W12_OP_xr.nc", engine = "netcdf4").rename({'__xarray_dataarray_variable__': 'tmpw'})
W31_OP_xr = xr.open_dataset(in_folder + "/3_Tower_new/Tower_xr/" + "W31_OP_xr.nc", engine = "netcdf4").rename({'__xarray_dataarray_variable__': 'tmpw'})

E06_RE_xr = xr.open_dataset(in_folder + "/3_Tower_new/Tower_xr/" + "E06_RE_xr.nc", engine = "netcdf4").rename({'__xarray_dataarray_variable__': 'tmpw'})
W12_RE_xr = xr.open_dataset(in_folder + "/3_Tower_new/Tower_xr/" + "W12_RE_xr.nc", engine = "netcdf4").rename({'__xarray_dataarray_variable__': 'tmpw'})
W31_RE_xr = xr.open_dataset(in_folder + "/3_Tower_new/Tower_xr/" + "W31_RE_xr.nc", engine = "netcdf4").rename({'__xarray_dataarray_variable__': 'tmpw'})

E06_EX_xr = xr.open_dataset(in_folder + "/3_Tower_new/Tower_xr/" + "E06_EX_xr.nc", engine = "netcdf4").rename({'__xarray_dataarray_variable__': 'tmpw'})
W12_EX_xr = xr.open_dataset(in_folder + "/3_Tower_new/Tower_xr/" + "W12_EX_xr.nc", engine = "netcdf4").rename({'__xarray_dataarray_variable__': 'tmpw'})
W31_EX_xr = xr.open_dataset(in_folder + "/3_Tower_new/Tower_xr/" + "W31_EX_xr.nc", engine = "netcdf4").rename({'__xarray_dataarray_variable__': 'tmpw'})

# RR_section = slice("2021-05-07T23:51:30", "2021-05-08T00:20:00")
RR_section = slice("2021-05-07T23:41:00", "2021-05-08T00:20:00")
# E06_OP_xr = E06_OP_xr.sel(time = W12_OP_xr.time)
E06_OP_xr = E06_OP_xr.sel(time = RR_section)
W12_OP_xr = W12_OP_xr.sel(time = RR_section)
W31_OP_xr = W31_OP_xr.sel(time = RR_section)

#%% interpolate the temp at 10.7m
def T_10dot5(E06_RE2, E06_H_prof):
    exp_T = np.mean(E06_RE2, axis=1)[15:36]
    exp_H = E06_H_prof[15:36]
    HH = np.arange(5, 12, 0.1)
    # plt.figure()
    s = InterpolatedUnivariateSpline(exp_H, exp_T, k = 1)
    y = s(HH)
#     plt.plot(y, HH)
    print("the temperature at 10.7 height is " + str(s(10.7)))
    return(s(10.7).item())

E06_T10 = T_10dot5(E06_RE_xr.tmpw, E06_RE_xr.z)
W12_T10 = T_10dot5(W12_RE_xr.tmpw, W12_RE_xr.z)
W31_T10 = T_10dot5(W31_RE_xr.tmpw, W31_RE_xr.z)

T10 = np.mean([E06_T10, W12_T10, W31_T10])

#%% first we start with TH plot
# basic settings
base_dir = os.path.dirname(os.path.realpath(__file__))
Case = base_dir +  "/W12"
# # Case_REF = base_dir +  ""
# sim_time = 3000
# TimeF = pd.read_pickle(base_dir + "/timeframe")
# Height = np.arange(0, 11)

#%% 
Ta_slice = np.load(Case + "/W12_Ta.npy")
#%%
TH = np.mean(Ta_slice, axis= 1)
# Tprof = np.mean(TH, axis= 0)
Height = np.arange(11)

#%% Time averaged profile with 75% quantile
EXP_RE_section = slice("2021-05-07T23:00:00", "2021-05-07T23:10:00")
W12_RE_exp = W12_EX_xr.sel(time = EXP_RE_section)
quartile1_re2, quartile3_re2 = np.percentile(W12_RE_exp.tmpw, [25, 75], axis=1)

fig2 = plt.figure(figsize=(6, 6))
ax = fig2.add_subplot(1,1,1)
h6 = ax.plot(np.mean(W12_RE_exp.tmpw, axis=1), W12_RE_exp.z, "-", label="on exp",linewidth = 4, ms = 6)   
# for i in range(TH.shape[0]):
for i in range(1, TH.shape[0], 4):
    h4 = ax.plot(TH[i,:], Height, "-o", label="off sim",linewidth = 1, ms = 2)  
# h4 = ax.plot(TH[-16,:], Height, "-", label="off sim",linewidth = 1, ms = 2)  
ax.fill_betweenx(W12_RE_exp.z, quartile1_re2, quartile3_re2, alpha=0.3, facecolor = "grey")
labels = np.arange(0,12,1)
plt.yticks(labels, labels)
plt.setp(ax, yticks=labels, xlim = (-5,11), ylim = (0, 9.5))
ax.set_ylabel('Height [m]',fontsize=18)
ax.set_xlabel("Temperature [$^\circ C$]",fontsize=18, usetex = False)
ax.set_ylim(0, 12)
plt.axvline(x=E06_T10, linestyle="--", linewidth = 1.5, color="k", alpha = 0.6)
# x_ticks = np.append(ax.get_xticks(), E06_T10)
# ax.set_xticks(x_ticks)
ax.text(x = E06_T10-0.8, y = 0.1, s = "9.6", fontsize = 18, c = "k")
# ax.set_xlim(-5, 11)
# legend = plt.legend(loc='best', frameon=False,fontsize = 18)
# plt.title("Averaged temp profile E06 [05-07]",fontsize=18)
plt.grid(linestyle = ":")
plt.tight_layout()
plt.savefig("b_prof.png")

# %%
