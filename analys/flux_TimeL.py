# %% ---------------------------------------------------------- #
#               -- quadrant analysis of flux # 
# % ----------------------------------------------------------- #
#%%
from tkinter import font
import xarray as xr
from glob import glob
import numpy as np
import matplotlib.pyplot as plt

plt.rc('font', family='serif')
plt.rc('xtick', labelsize='18')
plt.rc('ytick', labelsize='18')
plt.rc('text', usetex=True)
# %%
file_dir = "/home/dai/Documents/Reserach/Fruit/Krabbendijke_experiment_2021_0507/sonic_data/"
nc_file = glob(file_dir + "*.nc")

W31_file = xr.open_dataset(nc_file[1])
W31 = W31_file.to_dataframe()

W12_file = xr.open_dataset(nc_file[0])
W12 = W12_file.to_dataframe()

#%% plot the time series of wind
# W12_1s_OP = W12.rolling("1S").mean().loc['2021-05-07 23:40:00':'2021-05-08 01:00:00']
# W31_1s_OP = W31.rolling("1S").mean().loc['2021-05-07 23:40:00':'2021-05-08 01:00:00']
# W12_1s_RE2 = W12.rolling("1S").mean().loc['2021-05-07 21:41:00':'2021-05-07 23:00:00']
# W31_1s_RE2 = W31.rolling("1S").mean().loc['2021-05-07 21:41:00':'2021-05-07 23:00:00']

W12l = W12.loc['2021-05-07 21:41:00':'2021-05-08 00:20:00']
W31l = W31.loc['2021-05-07 21:41:00':'2021-05-08 00:20:00']

AP = 60
# the center can only be used with fixed number
W12_AT = W12.rolling(window = 600, center = True).mean().loc['2021-05-07 21:41:00':'2021-05-08 00:20:00']
W31_AT = W31.rolling(window = 600, center = True).mean().loc['2021-05-07 21:41:00':'2021-05-08 00:20:00']
# W12_AT = W12.loc['2021-05-07 21:40:00':'2021-05-08 00:21:00'].rolling(window = "%sS"%(str(AP))).mean()
# W31_AT = W31.loc['2021-05-07 21:40:00':'2021-05-08 00:21:00'].rolling(window = "%sS"%(str(AP))).mean()
# W12_AT = W12_AT.set_index(W12l.index)
# W31_AT = W31_AT.set_index(W31l.index)
# W12_1T_OP = W12.rolling("%sS"%(str(AP))).mean().loc['2021-05-07 23:40:00':'2021-05-08 01:00:00']
# W31_1T_OP = W31.rolling("%sS"%(str(AP))).mean().loc['2021-05-07 23:40:00':'2021-05-08 01:00:00']
# W12_1T_RE2 = W12.rolling("%sS"%(str(AP))).mean().loc['2021-05-07 21:41:00':'2021-05-07 23:00:00']
# W31_1T_RE2 = W31.rolling("%sS"%(str(AP))).mean().loc['2021-05-07 21:41:00':'2021-05-07 23:00:00']
#%% plot the detrend process -- whole series
fig1 = plt.figure(figsize=(10, 6))
ax = fig1.add_subplot(1,1,1)

h1 = plt.plot(W12l.index, W12l.u+3, 'k-', label='Raw',linewidth = 2)
h2 = plt.plot(W12_AT.index, W12_AT.u+3, 'r-', label='rolling',linewidth = 2)
h2 = plt.plot(W12_AT.index, W12l.u - W12_AT.u, 'b-', label='fluc',linewidth = 2)
legend = plt.legend(loc='best', frameon=False,fontsize = 18)
plt.grid(linestyle = ':')
plt.tight_layout()
#%% calculate the flux wT
W12_OP = W12.loc['2021-05-07 23:41:00':'2021-05-08 00:20:00']
W31_OP = W31.loc['2021-05-07 23:41:00':'2021-05-08 00:20:00']
W12_RE2 = W12.loc['2021-05-07 22:20:00':'2021-05-07 23:00:00']
W31_RE2 = W31.loc['2021-05-07 22:20:00':'2021-05-07 23:00:00']

W12_AT_OP = W12.rolling(window = 600, center = True).mean().loc['2021-05-07 23:41:00':'2021-05-08 00:20:00']
W31_AT_OP = W31.rolling(window = 600, center = True).mean().loc['2021-05-07 23:41:00':'2021-05-08 00:20:00']
W12_AT_RE2 =W12.rolling(window = 600, center = True).mean().loc['2021-05-07 22:20:00':'2021-05-07 23:00:00']
W31_AT_RE2 =W31.rolling(window = 600, center = True).mean().loc['2021-05-07 22:20:00':'2021-05-07 23:00:00']
#%%
W12_fluc_OP =  W12_OP - W12_AT_OP
W31_fluc_OP =  W31_OP - W31_AT_OP
W12_fluc_RE2 =  W12_RE2 - W12_AT_RE2
W31_fluc_RE2 =  W31_RE2 - W31_AT_RE2

W12_OP["wT"] = (W12_fluc_OP.w) * (W12_fluc_OP.Ts)
W12_OP["uw"] = (W12_fluc_OP.w) * (W12_fluc_OP.u)
W31_OP["wT"] = (W31_fluc_OP.w) * (W31_fluc_OP.Ts)
W31_OP["uw"] = (W31_fluc_OP.w) * (W31_fluc_OP.u)
W12_RE2["wT"] = (W12_fluc_RE2.w) * (W12_fluc_RE2.Ts)
W12_RE2["uw"] = (W12_fluc_RE2.w) * (W12_fluc_RE2.u)
W31_RE2["wT"] = (W31_fluc_RE2.w) * (W31_fluc_RE2.Ts)
W31_RE2["uw"] = (W31_fluc_RE2.w) * (W31_fluc_RE2.u)

#%%
fig1 = plt.figure(figsize=(6, 6))
ax = fig1.add_subplot(1,1,1)
# h1 = plt.plot(W12_OP.index, W12_OP["wT"], 'r--', label='wT',linewidth = 1)
# h1 = plt.plot(W12_OP.index, W12_OP["uw"], 'k--', label='wu',linewidth = 1)
h1 = plt.plot(W12_OP.index, W12_OP["M"], 'r--', label='u',linewidth = 1)
h1 = plt.plot(W12_OP.index, W12_OP["v"], 'k--', label='v',linewidth = 1)
# labels = np.arange(0,9.6,1)
# plt.yticks(labels, labels)
legend = plt.legend(loc='best', frameon=False,fontsize = 18)
plt.grid(linestyle = ':')
plt.tight_layout()
# plt.savefig('output_plots/...')
# %% phase average flux data
fl_t0 = 90
fl_dt = 2880

W12_flux_AV = W12_OP.iloc[fl_t0+fl_dt*0:fl_t0+fl_dt*1].values
for i in range(1, 5):
    W12_flux_AV = W12_flux_AV + W12_OP.iloc[fl_t0+fl_dt*i:fl_t0+fl_dt*(i+1)].values
W12_flux_AV = W12_flux_AV/5


#%% load the data from averaged tower
W12_OP_xr = xr.load_dataset("W12_OP_xr.nc", engine = "netcdf4").rename({'__xarray_dataarray_variable__': 'tmpw'}).tmpw
#%%
t0_W12 = 9
l_t_W12 = 288

W12_AV = W12_OP_xr.values[:,t0_W12+l_t_W12*0:t0_W12+l_t_W12*1]
for i in range(1, 5):
    W12_AV = W12_AV + W12_OP_xr.values[:,t0_W12+l_t_W12*i:t0_W12+l_t_W12*(i+1)]
W12_AV = W12_AV/5

# %%
W12_H_prof = np.linspace(0, 9.35, 37)
fig, (ax1) = plt.subplots(1, 1, figsize = (8, 6))
T22 = np.arange(-144, 144, 1)
TT = np.arange(-144, 144, 0.1)
# ax1.plot(TT, W12_flux_AV[:,7], "k+")
# ax1.plot(TT, W12_flux_AV[:,7]+3, "k-+", label = "wT")
# ax1.plot(TT, W12_flux_AV[:,8], "k-+", label = "wu")
ax1.plot(TT, pd.Series(W12_flux_AV[:,7]).rolling(20).mean(), "k-+", label = "wu")
# tp1 = ax1.contourf(T22, W12_H_prof, W12_AV, levels = np.arange(3, 9, 0.1), extend = "max", cmap = "gnuplot2")
# tpp1 = ax1.contour(T22, W12_H_prof, W12_AV, levels = np.arange(3, 9, 0.5))
# ax1.clabel(tpp1, fmt='%2.1f', colors='w', fontsize=14)
ax1.set_xlabel("Time relative to jet [s]", fontsize = 20)
plt.tight_layout()

# %%
