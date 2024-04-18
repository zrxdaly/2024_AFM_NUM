# %%
import glob
import os 
from tkinter import font
import xarray as xr
import glob
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import matplotlib.pyplot as plt
# plt.rc('font', family='serif')
plt.rc('xtick', labelsize='20')
plt.rc('ytick', labelsize='20')
plt.rc('text', usetex=True)

import scienceplots
plt.style.use(['science', 'grid', 'scatter'])
#%%
base_dir = os.path.dirname(os.path.realpath(__file__))
# %% U
Ufiles = sorted(glob.glob(base_dir + "/W12/W12_Ut=*x=528"))
USlice = np.fromfile(Ufiles[0], dtype = float).reshape(599, 11, 11)
for i in range(1, len(Ufiles)-1):
    print(Ufiles[i])
    USlice0 = np.fromfile(Ufiles[i], dtype = float).reshape(600, 11, 11)
    USlice = np.append(USlice, USlice0, axis = 0)
USlice0 = np.fromfile(Ufiles[len(Ufiles)-1], dtype = float).reshape(1, 11, 11)
USlice = np.append(USlice, USlice0, axis = 0)
# %% V
Vfiles = sorted(glob.glob(base_dir + "/W12/W12_Vt=*x=528"))
VSlice = np.fromfile(Vfiles[0], dtype = float).reshape(599, 11, 11)
for i in range(1, len(Vfiles)-1):
    print(Vfiles[i])
    VSlice0 = np.fromfile(Vfiles[i], dtype = float).reshape(600, 11, 11)
    VSlice = np.append(VSlice, VSlice0, axis = 0)
VSlice0 = np.fromfile(Vfiles[len(Vfiles)-1], dtype = float).reshape(1, 11, 11)
VSlice = np.append(VSlice, VSlice0, axis = 0)
# %% W
Wfiles = sorted(glob.glob(base_dir + "/W12/W12_Wt=*x=528"))
WSlice = np.fromfile(Wfiles[0], dtype = float).reshape(599, 11, 11)
for i in range(1, len(Wfiles)-1):
    print(Wfiles[i])
    WSlice0 = np.fromfile(Wfiles[i], dtype = float).reshape(600, 11, 11)
    WSlice = np.append(WSlice, WSlice0, axis = 0)
WSlice0 = np.fromfile(Wfiles[len(Wfiles)-1], dtype = float).reshape(1, 11, 11)
WSlice = np.append(WSlice, WSlice0, axis = 0)
#%%
U_h3 = (USlice[:,4,4] + USlice[:,4,4])/2.
V_h3 = (VSlice[:,4,4] + VSlice[:,4,4])/2.
W_h3 = (WSlice[:,4,4] + WSlice[:,4,4])/2.
M_h3 = np.sqrt(U_h3**2 + V_h3**2 + W_h3**2)
#%% phased average of the simulation velocity
M_h3_OP = M_h3[600:]
sp = 706-600
PAV_M_OP_h3 = M_h3_OP[sp+288*0:sp+288*(0+1)]
for j in range(1, 7):
    print(j)
    PAV_M_OP_h3 = PAV_M_OP_h3 + M_h3_OP[sp+288*j:sp+288*(j+1)]
PAV_M_sim_h3 = PAV_M_OP_h3/7
#%% phased average of the measurement velocity
sonic_dir = "/home/dai/Documents/Reserach/Fruit/Krabbendijke_experiment_2021_0507/sonic_data/"
nc_file = glob.glob(sonic_dir + "*.nc")

W12_file = xr.open_dataset(nc_file[1])
W12 = W12_file.to_dataframe()

W12_10s_EX = W12.rolling("10S").mean().loc['2021-05-07 22:20:00':'2021-05-08 00:20:00']
# %% phased average n -- first rotation time window: 23:41:16 --- 23:46:04 :: W12_10s_EX.index[48761:48761+288*10]
sttime = 48761
PAV_M_EX_h3 = W12_10s_EX["M"].values[sttime+288*10*0:sttime+288*10*(0+1)]
for j in range(1, 8):
    print(j)
    PAV_M_EX_h3 = PAV_M_EX_h3 + W12_10s_EX["M"].values[sttime+288*10*j:sttime+288*10*(j+1)]

PAV_M_obs_h3 = PAV_M_EX_h3/8

# %% compare the simualtion and measurement
tm = np.arange(0, 288, 0.1)
ts = np.arange(0, 288, 1)
fig1 = plt.figure(figsize=(7, 6))
ax = fig1.add_subplot(1,1,1)
ax.plot(ts, PAV_M_sim_h3, '-+', linewidth = 2, label = "Simulation")
ax.plot(tm, PAV_M_obs_h3, '-o', linewidth = 2, label = "Measurement")
ax.set_ylabel("M [$\mathrm{m~s^{-1}}$]",fontsize=18)
ax.set_xlabel("Time [s]",fontsize=18)
plt.grid(linestyle = ':')
plt.setp(ax, ylim = (0, 7), xlim = (0, 288))
legend = plt.legend(loc='upper left', frameon=False,fontsize = 18)
plt.tight_layout()
plt.savefig("W12/W12_M_OP.pdf")
# %% check some statistic 
M_h3_OP = M_h3[600:]
print(np.mean(M_h3_OP), np.std(M_h3_OP), np.max(M_h3_OP), np.min(M_h3_OP))

M_h3_RE = M_h3[:600]
print(np.mean(M_h3_RE), np.std(M_h3_RE), np.max(M_h3_RE), np.min(M_h3_RE))
# %%
W12_10s_RE = W12.rolling("10S").mean().loc['2021-05-07 22:50:00':'2021-05-07 23:00:00']
W12_10s_OP = W12.rolling("10S").mean().loc['2021-05-07 23:41:00':'2021-05-08 00:20:00']

print(np.mean(W12_10s_OP["M"]), np.std(W12_10s_OP["M"]), np.max(W12_10s_OP["M"]), np.min(W12_10s_OP["M"]))
print(np.mean(W12_10s_RE["M"]), np.std(W12_10s_RE["M"]), np.max(W12_10s_RE["M"]), np.min(W12_10s_RE["M"]))

# %%
