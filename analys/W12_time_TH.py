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
import glob 

plt.rc('font', family='serif')
plt.rc('xtick', labelsize='20')
plt.rc('ytick', labelsize='20')
plt.rc('text', usetex=True)
# import scienceplots
import scienceplots
plt.style.use(['science', 'grid', 'scatter'])

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
# Case_REF = base_dir +  ""
sim_time = 3000
TimeF = pd.read_pickle(base_dir + "/timeframe")
Height = np.arange(0, 11)

#%% 
Bfiles = sorted(glob.glob(base_dir + "/W12/W12_Bt=*x=528"))
BSlice = np.fromfile(Bfiles[0], dtype = float).reshape(599, 11, 11)
for i in range(1, len(Bfiles)-1):
    print(Bfiles[i])
    BSlice0 = np.fromfile(Bfiles[i], dtype = float).reshape(600, 11, 11)
    BSlice = np.append(BSlice, BSlice0, axis = 0)
BSlice0 = np.fromfile(Bfiles[len(Bfiles)-1], dtype = float).reshape(1, 11, 11)
BSlice = np.append(BSlice, BSlice0, axis = 0)
BSlice = BSlice/9.81*273.15
#%%
TH = BSlice[:,:,4]
T_shift = 130
fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize = (12, 8), sharex=True)
# ax3.plot(TimeF.index[T_shift:sim_time+T_shift], TH[:, 4], 'b-', label=r'$T_{1m}$',linewidth = 2)
ax3.plot(TimeF.index[T_shift:sim_time+T_shift], (TH[:, 1]+TH[:, 2])/2, 'b-', label=r'$T_{1m}$',linewidth = 2)
ax3.plot(TimeF.index[T_shift:sim_time+T_shift], TH[:, 9], 'r-', label=r'$T_{9m}$',linewidth = 2)
tp2 = ax2.contourf(TimeF.index[T_shift:sim_time+T_shift], Height, TH.T, levels = np.arange(0, 11, 0.1), cmap = "plasma")

tp1 = ax1.contourf(W12_EX_xr.time.values, W12_EX_xr.z.values, W12_EX_xr.tmpw.values, levels = np.arange(0, 10, 0.1), cmap = "plasma")
ax3.plot(W12_EX_xr.time.values, W12_EX_xr.tmpw[4,:], 'b*', label=r'$T_{1m}$',markersize = 4.5)
ax3.plot(W12_EX_xr.time.values, W12_EX_xr.tmpw[-2,:], 'r*', label=r'$T_{9m}$',markersize =4.5)

for c in tp1.collections:
    c.set_edgecolor("face")
for c in tp2.collections:
    c.set_edgecolor("face")

fig.subplots_adjust(bottom=0.1, top=0.9, left=0.1, right=0.8,
                wspace=0.02, hspace=0.05)
cb_ax = fig.add_axes([0.83, 0.1, 0.02, 0.8])
cbar = plt.colorbar(tp2, cax=cb_ax)
cbar.ax.set_ylabel('Air temperature [$\mathrm{^\circ C}$]',fontsize = 20)

ax3.set_ylabel("Temperature [$\mathrm{^\circ C}$]",fontsize=18)
ax2.set_ylabel("Height [m]",fontsize=18)
ax1.set_ylabel("Height [m]",fontsize=18)
ax3.set_xlabel("Time [hh:mm]",fontsize=18)
labels2 = np.arange(0,10,2)
labels1 = np.arange(0,12,2)
plt.setp(ax3, yticks = labels1, ylim = (0, 12), xlim = [pd.to_datetime("2021-05-07T23:34:00"), pd.to_datetime("2021-05-08T00:20:00")])
plt.setp((ax1, ax2), yticks = labels2, ylim = (0, 9.4), xlim = [pd.to_datetime("2021-05-07T23:34:00"), pd.to_datetime("2021-05-08T00:20:00")])

ax1.text(x = datetime.datetime(2021,5,8,0,18,10), y = 8.0, s = r"\textbf{a)}", fontsize = 20, c = "k")
ax2.text(x = datetime.datetime(2021,5,8,0,18,10), y = 8.0, s = r"\textbf{b)}", fontsize = 20, c = "k")
ax3.text(x = datetime.datetime(2021,5,8,0,18,10), y = 10.5, s = r"\textbf{c)}", fontsize = 20, c = "k")

legend = ax3.legend(loc='upper left', frameon=False,fontsize = 18, ncols = 5)
hh_mm = DateFormatter('%H:%M')
ax2.xaxis.set_major_formatter(hh_mm)
plt.savefig(Case+ "/Time_TH_exp_num.pdf", bbox_inches='tight')
# %% get the op and RE
TH_OP = BSlice[600:, :, 4]
TH_REF = BSlice[500:600, :, 4]
Height = np.arange(0, 11)
#%% Time averaged profile with 75% quantile
EXP_RE_section = slice("2021-05-07T22:50:00", "2021-05-07T23:00:00")
EXP_OP_section = slice("2021-05-07T23:41:00", "2021-05-08T00:20:00")
W12_RE_exp = W12_EX_xr.sel(time = EXP_RE_section)
W12_OP_exp = W12_EX_xr.sel(time = EXP_OP_section)

quartile1_op, quartile3_op = np.percentile(W12_OP_exp.tmpw, [25, 75], axis=1)
quartile1_re2, quartile3_re2 = np.percentile(W12_RE_exp.tmpw, [25, 75], axis=1)

quartile1_TH_op, quartile3_TH_op = np.percentile(TH_OP, [25, 75], axis=0)
quartile1_TH_ref, quartile3_TH_ref = np.percentile(TH_REF, [10, 90], axis=0)

fig2 = plt.figure(figsize=(6, 6))
ax = fig2.add_subplot(1,1,1)
h5 = ax.plot(np.mean(W12_OP_exp.tmpw, axis=1), W12_OP_exp.z, "r+", label="on exp",linewidth = 2, ms = 6)  
h6 = ax.plot(np.mean(W12_RE_exp.tmpw, axis=1), W12_RE_exp.z, "k+", label="off exp",linewidth = 2, ms = 6)  
h3 = ax.plot(np.mean(TH_OP, axis=0), Height, "r-o", label="on sim",linewidth = 2, ms = 6)  
h4 = ax.plot(np.mean(TH_REF, axis=0), Height, "k-o", label="off sim",linewidth = 2, ms = 6)  

ax.fill_betweenx(W12_OP_exp.z, quartile1_op, quartile3_op, alpha=0.3, facecolor = "r")
ax.fill_betweenx(W12_RE_exp.z, quartile1_re2, quartile3_re2, alpha=0.3, facecolor = "grey")

ax.fill_betweenx(Height, quartile1_TH_op, quartile3_TH_op, alpha=0.3, facecolor = "r")
ax.fill_betweenx(Height, quartile1_TH_ref, quartile3_TH_ref, alpha=0.3, facecolor = "grey")

labels = np.arange(0,9.6,1)
plt.yticks(labels, labels)
plt.setp(ax, yticks=labels, xlim = (-1,11), ylim = (0, 9.5))
ax.set_ylabel('Height [m]',fontsize=18)
ax.set_xlabel("Air temperature [$\mathrm{^\circ C}$]",fontsize=18)
ax.set_ylim(0, 9.5)
plt.axvline(x=E06_T10, linestyle="--", linewidth = 1.5, color="k", alpha = 0.6)
# x_ticks = np.append(ax.get_xticks(), E06_T10)
# ax.set_xticks(x_ticks)
ax.text(x = E06_T10-0.8, y = 0.1, s = "9.6", fontsize = 18, c = "k")
ax.set_xlim(-1, 11)
legend = plt.legend(loc='best', frameon=False,fontsize = 18)
# plt.title("Averaged temp profile E06 [05-07]",fontsize=18)
plt.grid(linestyle = ":")
plt.tight_layout()
plt.savefig(Case + "/T_profile_exp.pdf")

#%% # ------------adding the wind time series to the third plot %% U--------------------------- #
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
U_h3 = USlice[:,4,4]
V_h3 = VSlice[:,4,4]
W_h3 = WSlice[:,4,4]
M_h3 = np.sqrt(U_h3**2 + V_h3**2 + W_h3**2)
#%%
sonic_dir = "/home/dai/Documents/Reserach/Fruit/Krabbendijke_experiment_2021_0507/sonic_data/"
nc_file = glob.glob(sonic_dir + "*.nc")

W12_file = xr.open_dataset(nc_file[1])
W12 = W12_file.to_dataframe()

W12_EX = W12.loc['2021-05-07 22:20:00':'2021-05-08 00:20:00']
#%%
# sim_time = 3000
# T_shift = 10
TH = BSlice[:,:,4]
T_shift = 150
fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, figsize = (12, 9), sharex=True)
# ax3.plot(TimeF.index[T_shift:sim_time+T_shift], TH[:, 1], 'b-', label=r'$T_{1m}$',linewidth = 2)
ax3.plot(TimeF.index[T_shift:sim_time+T_shift], (TH[:, 1]+TH[:, 2])/2, 'b-', label=r'$T_{1m}$',linewidth = 2)
ax3.plot(TimeF.index[T_shift:sim_time+T_shift], TH[:, 9], 'r-', label=r'$T_{9m}$',linewidth = 2)
tp2 = ax2.contourf(TimeF.index[T_shift:sim_time+T_shift], Height, TH.T, levels = np.arange(0, 11, 0.1), cmap = "plasma")

tp1 = ax1.contourf(W12_EX_xr.time.values, W12_EX_xr.z.values, W12_EX_xr.tmpw.values, levels = np.arange(0, 10, 0.1), cmap = "plasma")
ax3.plot(W12_EX_xr.time.values, W12_EX_xr.tmpw[4,:], 'b*', label=r'$T_{1m}$',markersize = 4.5)
ax3.plot(W12_EX_xr.time.values, W12_EX_xr.tmpw[-2,:], 'r*', label=r'$T_{9m}$',markersize =4.5)

ax4.plot(W12_EX.index, W12_EX.M, '*', label='M obs',markersize =1)
ax4.plot(TimeF.index[T_shift:sim_time+T_shift], M_h3, '+', label='M sim',markersize =2)

for c in tp1.collections:
    c.set_edgecolor("face")
for c in tp2.collections:
    c.set_edgecolor("face")

fig.subplots_adjust(bottom=0.1, top=0.9, left=0.1, right=0.8,
                wspace=0.02, hspace=0.05)
cb_ax = fig.add_axes([0.83, 0.1, 0.02, 0.8])
cbar = plt.colorbar(tp2, cax=cb_ax)
cbar.ax.set_ylabel('Air temperature [$\mathrm{^\circ C}$]',fontsize = 20)

ax3.set_ylabel("Temperature [$\mathrm{^\circ C}$]",fontsize=18)
ax4.set_ylabel("M [$\mathrm{m~s^{-1}}$]",fontsize=18)
ax2.set_ylabel("Height [m]",fontsize=18)
ax1.set_ylabel("Height [m]",fontsize=18)
ax4.set_xlabel("Time [hh:mm]",fontsize=18)
labels2 = np.arange(0,10,2)
labels1 = np.arange(0,12,2)
labelsU = np.arange(-4,15,4)
plt.setp(ax3, yticks = labels1, ylim = (0, 12), xlim = [pd.to_datetime("2021-05-07T23:34:00"), pd.to_datetime("2021-05-08T00:20:00")])
plt.setp(ax4, yticks = labelsU, ylim = (-4, 15), xlim = [pd.to_datetime("2021-05-07T23:34:00"), pd.to_datetime("2021-05-08T00:20:00")])
plt.setp((ax1, ax2), yticks = labels2, ylim = (0, 9.4), xlim = [pd.to_datetime("2021-05-07T23:34:00"), pd.to_datetime("2021-05-08T00:20:00")])

ax1.text(x = datetime.datetime(2021,5,8,0,18,10), y = 8.0, s = r"\textbf{a)}", fontsize = 20, c = "k")
ax2.text(x = datetime.datetime(2021,5,8,0,18,10), y = 8.0, s = r"\textbf{b)}", fontsize = 20, c = "k")
ax3.text(x = datetime.datetime(2021,5,8,0,18,10), y = 10.5, s = r"\textbf{c)}", fontsize = 20, c = "k")

legend = ax3.legend(loc='upper left', frameon=False,fontsize = 18, ncols = 5)
hh_mm = DateFormatter('%H:%M')
ax2.xaxis.set_major_formatter(hh_mm)
# plt.savefig(Case+ "/Time_TH_exp_num.pdf", bbox_inches='tight')
# %%
