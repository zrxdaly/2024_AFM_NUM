#%% this file plot the Ta and Tv
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from matplotlib.dates import DateFormatter
import datetime
import glob

# plt.rc('font', family='serif')
plt.rc('xtick', labelsize='20')
plt.rc('ytick', labelsize='20')
plt.rc('text', usetex=True)
import scienceplots
plt.style.use(['science', 'grid', 'scatter'])
# %%
Tva_file = pd.read_csv("20210507_Midden.csv", header = 1, names = ["i", "date", "Va", "Vv", "Vs", "Vt"], delimiter=",", parse_dates = [1], infer_datetime_format = True)
# %%
# Tva_file.set_index["date"]
Tva_file.set_index("date", inplace = True, drop = True)
# %% the leave air temperature from measurement
def V_to_C(V):
        a = 1.13075635e-3
        b = 2.33896902e-4
        c = 8.82996895e-8
        Rref = 20000
        Vex = 2.5
        Rntc = (Vex - V) / V * Rref
        T = 1 / (a + b * np.log(Rntc) + c * np.log(Rntc)**3) - 273.15 
        return T

#%%
Tva_file["Ta"] = V_to_C(Tva_file["Va"])
Tva_file["Tv"] = V_to_C(Tva_file["Vv"])

# %%
INT_time = slice("2021-05-07T23:40:00", "2021-05-08T00:20:00")
OP_T_section = slice("2021-05-07T23:40:00", "2021-05-08T00:20:00")
RE_T_section = slice("2021-05-07T22:20:00", "2021-05-07T23:00:00")

Tva_file_INT  = Tva_file.loc['2021-05-07 22:20':'2021-05-08 00:20']
Tva_file_OP  = Tva_file.loc['2021-05-07 23:27:30':'2021-05-08 00:20:00']
Tva_file_RE  = Tva_file.loc['2021-05-07 22:20':'2021-05-07 22:40']

# %% load the simulation data
# basic settings
base_dir = os.path.dirname(os.path.realpath(__file__))
Case = base_dir +  "/W12"
TimeF = pd.read_pickle(base_dir + "/timeframe")

#%% 
timeV = 40
Bfiles = sorted(glob.glob(base_dir + "/W12/W12_Bt=*x=532"))
BSlice = np.fromfile(Bfiles[0], dtype = float).reshape(599, 11, 11)
for i in range(1, len(Bfiles)-1):
    print(Bfiles[i])
    BSlice0 = np.fromfile(Bfiles[i], dtype = float).reshape(600, 11, 11)
    BSlice = np.append(BSlice, BSlice0, axis = 0)
BSlice0 = np.fromfile(Bfiles[len(Bfiles)-1], dtype = float).reshape(1, 11, 11)
BSlice = np.append(BSlice, BSlice0, axis = 0)
BSlice = BSlice/9.81*273.15
#%%
TVfiles = sorted(glob.glob(base_dir + "/W12/W12_TVt=*x=532"))
TVSlice = np.fromfile(TVfiles[0], dtype = float).reshape(599, 11, 11)
for i in range(1, len(TVfiles)-1):
    print(TVfiles[i])
    TVSlice0 = np.fromfile(TVfiles[i], dtype = float).reshape(600, 11, 11)
    TVSlice = np.append(TVSlice, TVSlice0, axis = 0)
TVSlice0 = np.fromfile(TVfiles[len(TVfiles)-1], dtype = float).reshape(1, 11, 11)
TVSlice = np.append(TVSlice, TVSlice0, axis = 0)
TVSlice = TVSlice-273.15
#%%
SIM_time = 3000
# hei = 2
ind = 3
# TV_y15 = (TVSlice[:,2,ind] + TVSlice[:,1,ind])/2
# B_y15  = (BSlice[:,2,ind] + BSlice[:,1,ind])/2
TV_y15 = TVSlice[:,2,ind]
B_y15  = BSlice[:,2,ind] 
# TV_y15 = TVSlice[:,3,ind]
# B_y15  = BSlice[:,3,ind] 
# TV_y15 = TVSlice[:,1,ind]
# B_y15  = BSlice[:,1,ind] 
#%%
fig, (ax1, ax2) = plt.subplots(2, 1, figsize = (10, 6), sharex=True)

tp1 = ax1.plot(Tva_file_OP.index[:SIM_time], Tva_file_OP["Ta"][:SIM_time], 'k-', label=r'$\mathrm{T_{a}}$', linewidth = 4)
tp2 = ax1.plot(Tva_file_OP.index[:SIM_time], Tva_file_OP["Tv"][:SIM_time], 'g-', label=r'$\mathrm{T_{v}}$', linewidth = 4)
tp3 = ax2.plot(Tva_file_OP.index[timeV:SIM_time+timeV], B_y15, 'k-', label=r'$\mathrm{T_{a}}$', linewidth = 4)
tp4 = ax2.plot(Tva_file_OP.index[timeV:SIM_time+timeV], TV_y15, 'g-', label=r'$\mathrm{T_{p}}$', linewidth = 4)
ax1.set_ylabel("Temperature [$\mathrm{^\circ C}$]",fontsize=18)
ax2.set_ylabel("Temperature [$\mathrm{^\circ C}$]",fontsize=18)
ax2.set_xlabel("Time [hh:mm]",fontsize=18)
labels2 = np.arange(0, 10, 2)
labels1 = np.arange(0, 10, 2)
plt.setp(ax1, yticks = labels2, ylim = (0, 10))
plt.setp(ax2, yticks = labels1, ylim = (0, 10), xlim = [pd.to_datetime("2021-05-07T23:30:00"), pd.to_datetime("2021-05-08T00:17:30")])
legend = plt.legend(loc='best', frameon=False,fontsize = 18)
hh_mm = DateFormatter('%H:%M')
ax2.xaxis.set_major_formatter(hh_mm)
ax1.text(x = datetime.datetime(2021,5,8,00,7), y = 1, s = r"\textbf{a)} measurement", fontsize = 20, c = "k")
ax2.text(x = datetime.datetime(2021,5,8,00,8), y = 1, s = r"\textbf{b)} simulation", fontsize = 20, c = "k")
# ax1.grid(linestyle = ':')
# ax2.grid(linestyle = ':')
plt.tight_layout()
plt.savefig("W12/TV_TA.pdf")

# %%
