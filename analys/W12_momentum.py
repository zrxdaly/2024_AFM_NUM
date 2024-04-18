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

# %% # --------------------------------------------------------- #
file_dir = os.path.dirname(os.path.realpath(__file__))
Bfiles = sorted(glob.glob(file_dir + "/W12/W12_Bt=*x=528"))
BSlice = np.fromfile(Bfiles[0], dtype = float).reshape(599, 11, 11)
for i in range(1, len(Bfiles)-1):
    print(Bfiles[i])
    BSlice0 = np.fromfile(Bfiles[i], dtype = float).reshape(600, 11, 11)
    BSlice = np.append(BSlice, BSlice0, axis = 0)
BSlice0 = np.fromfile(Bfiles[len(Bfiles)-1], dtype = float).reshape(1, 11, 11)
BSlice = np.append(BSlice, BSlice0, axis = 0)
BSlice = BSlice/9.81*273.15
# %% U
Ufiles = sorted(glob.glob(file_dir + "/W12/W12_Ut=*x=528"))
USlice = np.fromfile(Ufiles[0], dtype = float).reshape(599, 11, 11)
for i in range(1, len(Ufiles)-1):
    print(Ufiles[i])
    USlice0 = np.fromfile(Ufiles[i], dtype = float).reshape(600, 11, 11)
    USlice = np.append(USlice, USlice0, axis = 0)
USlice0 = np.fromfile(Ufiles[len(Ufiles)-1], dtype = float).reshape(1, 11, 11)
USlice = np.append(USlice, USlice0, axis = 0)
# %% V
Vfiles = sorted(glob.glob(file_dir + "/W12/W12_Vt=*x=528"))
VSlice = np.fromfile(Vfiles[0], dtype = float).reshape(599, 11, 11)
for i in range(1, len(Vfiles)-1):
    print(Vfiles[i])
    VSlice0 = np.fromfile(Vfiles[i], dtype = float).reshape(600, 11, 11)
    VSlice = np.append(VSlice, VSlice0, axis = 0)
VSlice0 = np.fromfile(Vfiles[len(Vfiles)-1], dtype = float).reshape(1, 11, 11)
VSlice = np.append(VSlice, VSlice0, axis = 0)
# %% W
Wfiles = sorted(glob.glob(file_dir + "/W12/W12_Wt=*x=528"))
WSlice = np.fromfile(Wfiles[0], dtype = float).reshape(599, 11, 11)
for i in range(1, len(Wfiles)-1):
    print(Wfiles[i])
    WSlice0 = np.fromfile(Wfiles[i], dtype = float).reshape(600, 11, 11)
    WSlice = np.append(WSlice, WSlice0, axis = 0)
WSlice0 = np.fromfile(Wfiles[len(Wfiles)-1], dtype = float).reshape(1, 11, 11)
WSlice = np.append(WSlice, WSlice0, axis = 0)
#%% plot the U V W T at certain height 
fig1 = plt.figure(figsize=(10, 4))
ax = fig1.add_subplot(1,1,1)
idx = 4
plt.plot(np.arange(3000), USlice[:,idx,4], 'r-', linewidth = 2)
plt.plot(np.arange(3000), VSlice[:,idx,4], 'b-', linewidth = 2)
plt.plot(np.arange(3000), WSlice[:,idx,4], 'g-', linewidth = 2)
plt.plot(np.arange(3000), BSlice[:,idx,4], 'k-', linewidth = 2)
# plt.xlim([-1,10])
# plt.ylim([0,600])
plt.grid(linestyle = ':')
plt.tight_layout()
#%% calculate simulated flux
U_h3 = USlice[:,4,4]
V_h3 = VSlice[:,4,4]
W_h3 = WSlice[:,4,4]
B_h3 = BSlice[:,4,4]

#%% taking average to calcuate the flux
# U_h3_RE = U_h3[:600]
# U_h3_OP = U_h3[600:]
# V_h3_RE = V_h3[:600]
# V_h3_OP = V_h3[600:]
# W_h3_RE = W_h3[:600]
# W_h3_OP = W_h3[600:]
# B_h3_RE = B_h3[:600]
# B_h3_OP = B_h3[600:]

# flc_U_RE = U_h3_RE - np.mean(U_h3_RE)
# flc_U_OP = U_h3_OP - np.mean(U_h3_OP)
# flc_V_RE = V_h3_RE - np.mean(V_h3_RE)
# flc_V_OP = V_h3_OP - np.mean(V_h3_OP)
# flc_W_RE = W_h3_RE - np.mean(W_h3_RE)
# flc_W_OP = W_h3_OP - np.mean(W_h3_OP)
# flc_B_RE = B_h3_RE - np.mean(B_h3_RE)
# flc_B_OP = B_h3_OP - np.mean(B_h3_OP)

# flx_UW_RE = flc_U_RE * flc_W_RE
# flx_UW_OP = flc_U_OP * flc_W_OP
# flx_VW_RE = flc_V_RE * flc_W_RE
# flx_VW_OP = flc_V_OP * flc_W_OP
# flx_BW_RE = flc_B_RE * flc_W_RE
# flx_BW_OP = flc_B_OP * flc_W_OP

#%% using moving average to calculate the flux
def moving_average(x, w):
    return np.convolve(x, np.ones(w), 'same') / w
#%%
flc_U = U_h3 - moving_average(U_h3, 60)
flc_V = V_h3 - moving_average(V_h3, 60)
flc_W = W_h3 - moving_average(W_h3, 60)
flc_B = B_h3 - moving_average(B_h3, 60)

flx_BW = flc_B * flc_W
flx_UW = flc_U * flc_W
flx_VW = flc_V * flc_W
 
# %%
sonic_dir = "/home/dai/Documents/Reserach/Fruit/Krabbendijke_experiment_2021_0507/sonic_data/"
nc_file = glob.glob(sonic_dir + "*.nc")

W31_file = xr.open_dataset(nc_file[1])
W31 = W31_file.to_dataframe()

W12_file = xr.open_dataset(nc_file[0])
W12 = W12_file.to_dataframe()

#%% plot the time series of wind
# W12l = W12.loc['2021-05-07 21:41:00':'2021-05-08 00:20:00']
# W31l = W31.loc['2021-05-07 21:41:00':'2021-05-08 00:20:00']

# AP = 60
# # rolling with 60S
# W12_AT = W12.rolling(window = 600, center = True).mean().loc['2021-05-07 22:20:00':'2021-05-08 00:20:00']
# W31_AT = W31.rolling(window = 600, center = True).mean().loc['2021-05-07 22:20:00':'2021-05-08 00:20:00']
# #%% plot the detrend process -- whole series
# fig1 = plt.figure(figsize=(12, 6))
# ax = fig1.add_subplot(1,1,1)

# h1 = plt.plot(W31l.index, W31l.u+3, 'k-', label='Raw',linewidth = 2)
# h2 = plt.plot(W31_AT.index, W31_AT.u+3, 'r-', label='rolling',linewidth = 2)
# h2 = plt.plot(W31_AT.index, W31l.u - W31_AT.u, 'b-', label='fluc',linewidth = 2)

# legend = plt.legend(loc='best', frameon=False,fontsize = 18)
# plt.grid(linestyle = ':')
# plt.tight_layout()
#%% calculate the flux wT
W12_EX = W12.loc['2021-05-07 22:20:00':'2021-05-08 00:20:00']
W31_EX = W31.loc['2021-05-07 22:20:00':'2021-05-08 00:20:00']
W12_AT_EX =W12.rolling(window = 600, center = True).mean().loc['2021-05-07 22:20:00':'2021-05-08 00:20:00']
W31_AT_EX =W31.rolling(window = 600, center = True).mean().loc['2021-05-07 22:20:00':'2021-05-08 00:20:00']

W12_fluc_EX =  W12_EX - W12_AT_EX
W31_fluc_EX =  W31_EX - W31_AT_EX

W12_fluc_EX["wT"] = (W12_fluc_EX.w) * (W12_fluc_EX.Ts)
W12_fluc_EX["uw"] = (W12_fluc_EX.w) * (W12_fluc_EX.u)

W31_fluc_EX["wT"] = (W31_fluc_EX.w) * (W31_fluc_EX.Ts)
W31_fluc_EX["uw"] = (W31_fluc_EX.w) * (W31_fluc_EX.u)
# %%
TimeF = pd.read_pickle(file_dir + "/timeframe")
#%%
sim_time = 3000
T_shift = 10
fig1 = plt.figure(figsize=(12, 6))
ax = fig1.add_subplot(1,1,1)
ax.plot(W12_fluc_EX.index, W12_fluc_EX.wT, 'r-', label='W12_OP')
plt.plot(TimeF.index[T_shift:sim_time+T_shift], flx_BW, 'b-', label='sim OP')
labels1 = np.arange(-1.5,1.5,0.5)
plt.setp(ax, yticks = labels1, ylim = (-1.5, 1.5), xlim = [pd.to_datetime("2021-05-07T23:33:00"), pd.to_datetime("2021-05-08T00:17:30")])
plt.grid(linestyle = ':')
plt.legend()
plt.tight_layout()
#%%
sim_time = 3000
T_shift = 10
fig1 = plt.figure(figsize=(12, 6))
ax = fig1.add_subplot(1,1,1)
ax.plot(W12_fluc_EX.index, W12_fluc_EX.uw, 'r-', label='W12_OP')
plt.plot(TimeF.index[T_shift:sim_time+T_shift], flx_UW, 'b-', label='sim OP')
labels1 = np.arange(-1.5,1.5,0.5)
plt.setp(ax, yticks = labels1, ylim = (-1.5, 1.5), xlim = [pd.to_datetime("2021-05-07T23:33:00"), pd.to_datetime("2021-05-08T00:17:30")])
plt.grid(linestyle = ':')
plt.legend()
plt.tight_layout()
# %%
