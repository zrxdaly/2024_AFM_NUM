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
# rolling with 60S
W12_AT = W12.rolling(window = 600, center = True).mean().loc['2021-05-07 21:41:00':'2021-05-08 00:20:00']
W31_AT = W31.rolling(window = 600, center = True).mean().loc['2021-05-07 21:41:00':'2021-05-08 00:20:00']

# W12_1T_OP = W12.rolling("%sS"%(str(AP))).mean().loc['2021-05-07 23:40:00':'2021-05-08 01:00:00']
# W31_1T_OP = W31.rolling("%sS"%(str(AP))).mean().loc['2021-05-07 23:40:00':'2021-05-08 01:00:00']
# W12_1T_RE2 = W12.rolling("%sS"%(str(AP))).mean().loc['2021-05-07 21:41:00':'2021-05-07 23:00:00']
# W31_1T_RE2 = W31.rolling("%sS"%(str(AP))).mean().loc['2021-05-07 21:41:00':'2021-05-07 23:00:00']
#%% plot the detrend process -- whole series
fig1 = plt.figure(figsize=(6, 6))
ax = fig1.add_subplot(1,1,1)

h1 = plt.plot(W12l.index, W12l.u+3, 'k-', label='Raw',linewidth = 2)
h2 = plt.plot(W12_AT.index, W12_AT.u+3, 'r-', label='rolling',linewidth = 2)
h2 = plt.plot(W12_AT.index, W12l.u - W12_AT.u, 'b-', label='fluc',linewidth = 2)

legend = plt.legend(loc='best', frameon=False,fontsize = 18)
plt.grid(linestyle = ':')
plt.tight_layout()

#%% plot the detrend process -- whole series
fig1 = plt.figure(figsize=(6, 6))
ax = fig1.add_subplot(1,1,1)

h1 = plt.plot(W31l.index, W31l.u+3, 'k-', label='Raw',linewidth = 2)
h2 = plt.plot(W31_AT.index, W31_AT.u+3, 'r-', label='rolling',linewidth = 2)
h2 = plt.plot(W31_AT.index, W31l.u - W31_AT.u, 'b-', label='fluc',linewidth = 2)

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

W12_fluc_OP =  W12_OP - W12_AT_OP
W31_fluc_OP =  W31_OP - W31_AT_OP
W12_fluc_RE2 =  W12_RE2 - W12_AT_RE2
W31_fluc_RE2 =  W31_RE2 - W31_AT_RE2

W12_fluc_OP["wT"] = (W12_fluc_OP.w) * (W12_fluc_OP.Ts)
W12_fluc_OP["uw"] = (W12_fluc_OP.w) * (W12_fluc_OP.u)

W31_fluc_OP["wT"] = (W31_fluc_OP.w) * (W31_fluc_OP.Ts)
W31_fluc_OP["uw"] = (W31_fluc_OP.w) * (W31_fluc_OP.u)

W12_fluc_RE2["wT"] = (W12_fluc_RE2.w) * (W12_fluc_RE2.Ts)
W12_fluc_RE2["uw"] = (W12_fluc_RE2.w) * (W12_fluc_RE2.u)

W31_fluc_RE2["wT"] = (W31_fluc_RE2.w) * (W31_fluc_RE2.Ts)
W31_fluc_RE2["uw"] = (W31_fluc_RE2.w) * (W31_fluc_RE2.u)

# %% ---------------------------------------------------------- #
# -- got the flux variable from 3D sonic, now put them into different hole # 
# % ----------------------------------------------------------- #
#%%  first show the correlation
print("The correlation between w and T for W12 duing OP", np.corrcoef(W12_fluc_OP.w, W12_fluc_OP.Ts)[0,1])
print("The correlation between w and T for W31 duing OP", np.corrcoef(W31_fluc_OP.w, W31_fluc_OP.Ts)[0,1])
print("The correlation between w and T for W12 duing RE2", np.corrcoef(W12_fluc_RE2.w, W12_fluc_RE2.Ts)[0,1])
print("The correlation between w and T for W31 duing RE2", np.corrcoef(W31_fluc_RE2.w, W31_fluc_RE2.Ts)[0,1])
#%%
print("The correlation between w and u for W12 duing OP", np.corrcoef(W12_fluc_OP.w, W12_fluc_OP.u)[0,1])
print("The correlation between w and u for W31 duing OP", np.corrcoef(W31_fluc_OP.w, W31_fluc_OP.u)[0,1])
print("The correlation between w and u for W12 duing RE2", np.corrcoef(W12_fluc_RE2.w, W12_fluc_RE2.u)[0,1])
print("The correlation between w and u for W31 duing RE2", np.corrcoef(W31_fluc_RE2.w, W31_fluc_RE2.u)[0,1])
#%%  try to define the flux into four zone
def QH_analysis(Fluc_OP, H):
    # normalized with std()
    thres = np.std(Fluc_OP.w) * np.std(Fluc_OP.Ts)
    # thres_H = H * thres
    # normalized with mean w'T'
    thres_H = H * np.mean(np.abs(Fluc_OP['wT']))
    # SUM_FLUX = np.mean([np.abs(W12_fluc_OP['wT']).sum(), np.abs(W12_fluc_RE2['wT']).sum(), np.abs(W31_fluc_OP['wT']).sum(), np.abs(W31_fluc_RE2['wT']).sum()])
    # w at y and T at x
    # w>0, T>0
    Q1_W12 = Fluc_OP[(Fluc_OP["w"]>0)&(Fluc_OP["Ts"]>0)&(np.abs(Fluc_OP['wT'])>=thres_H)]
    # w<0, T>0
    Q4_W12 = Fluc_OP[(Fluc_OP["w"]<0)&(Fluc_OP["Ts"]>0)&(np.abs(Fluc_OP['wT'])>=thres_H)]
    # w<0, T<0
    Q3_W12 = Fluc_OP[(Fluc_OP["w"]<0)&(Fluc_OP["Ts"]<0)&(np.abs(Fluc_OP['wT'])>=thres_H)]
    # w>0, T<0
    Q2_W12 = Fluc_OP[(Fluc_OP["w"]>0)&(Fluc_OP["Ts"]<0)&(np.abs(Fluc_OP['wT'])>=thres_H)]

    Q1_tiH = Q1_W12.shape[0] / Fluc_OP.shape[0]
    Q2_tiH = Q2_W12.shape[0] / Fluc_OP.shape[0]
    Q3_tiH = Q3_W12.shape[0] / Fluc_OP.shape[0]
    Q4_tiH = Q4_W12.shape[0] / Fluc_OP.shape[0]

    # normalized with std()
    # Q1_SiH = Q1_W12["wT"].sum() * Q1_tiH / thres
    # Q2_SiH = Q2_W12["wT"].sum() * Q2_tiH / thres
    # Q3_SiH = Q3_W12["wT"].sum() * Q3_tiH / thres
    # Q4_SiH = Q4_W12["wT"].sum() * Q4_tiH / thres

    # SiH in this case are with signs 
    Q1_SiH = Q1_W12["wT"].sum() / np.abs(Fluc_OP['wT']).sum()
    Q2_SiH = Q2_W12["wT"].sum() / np.abs(Fluc_OP['wT']).sum()
    Q3_SiH = Q3_W12["wT"].sum() / np.abs(Fluc_OP['wT']).sum()
    Q4_SiH = Q4_W12["wT"].sum() / np.abs(Fluc_OP['wT']).sum()
    return([Q1_tiH, Q2_tiH, Q3_tiH, Q4_tiH, Q1_SiH, Q2_SiH, Q3_SiH, Q4_SiH])
#%%
HP = np.arange(0, 20, 0.1)
Q_W12_OP = []
for h in HP:
    Q = QH_analysis(W12_fluc_OP, h)
    Q_W12_OP = np.append(Q_W12_OP, Q)
Q_W12_OP = Q_W12_OP.reshape((HP.shape[0], 8))

Q_W12_RE2 = []
for h in HP:
    Q = QH_analysis(W12_fluc_RE2, h)
    Q_W12_RE2 = np.append(Q_W12_RE2, Q)
Q_W12_RE2 = Q_W12_RE2.reshape((HP.shape[0], 8))


#%% quadrant plot of W12 OP
fig1 = plt.figure(figsize=(6, 6))
ax = fig1.add_subplot(1,1,1)

h0 = ax.plot( HP,  Q_W12_RE2[:, 4], 'k-+', label='Q1',linewidth = 2, ms = 4)
h1 = ax.plot(-HP, -Q_W12_RE2[:, 5], 'k-+', label='Q2',linewidth = 2, ms = 4)
h2 = ax.plot(-HP, -Q_W12_RE2[:, 6], 'k-+', label='Q3',linewidth = 2, ms = 4)
h3 = ax.plot( HP,  Q_W12_RE2[:, 7], 'k-+', label='Q4',linewidth = 2, ms = 4)

h0 = ax.plot( HP,  Q_W12_OP[:, 4], 'r-+', label='Q1',linewidth = 2, ms = 4)
h1 = ax.plot(-HP, -Q_W12_OP[:, 5], 'r-+', label='Q2',linewidth = 2, ms = 4)
h2 = ax.plot(-HP, -Q_W12_OP[:, 6], 'r-+', label='Q3',linewidth = 2, ms = 4)
h3 = ax.plot( HP,  Q_W12_OP[:, 7], 'r-+', label='Q4',linewidth = 2, ms = 4)
ax.set_title("W12 heat flux fraction", fontsize = 18)
# plt.yscale('symlog', linthresh=0.05)
# plt.xscale('symlog')
# legend = plt.legend(loc='best', frameon=False,fontsize = 18)
# plt.grid(linestyle = ':')
plt.tight_layout()
# plt.savefig('output_plots/...')

#%% example plots of QH analysis
fig1 = plt.figure(figsize=(6, 6))
ax = fig1.add_subplot(1,1,1)

plt.scatter(W12_fluc_OP.u, W12_fluc_OP.w, marker = "*")
xx = np.arange(0.01, 20, 0.1)
np.mean(np.abs(W12_fluc_OP['uw']))
# yy = 6 / xx
# yy2 = 3 / xx
yy3 = 0.3 / xx
# plt.plot(xx, -yy , 'r--', label=r'$\left|u^{\prime} w^{\prime}\right| = 20 \overline{\left|u^{\prime} w^{\prime}\right|}$',linewidth = 2)
# plt.plot(xx, -yy2, 'k--', label=r'$\left|u^{\prime} w^{\prime}\right| = 10 \overline{\left|u^{\prime} w^{\prime}\right|}$',linewidth = 2)
plt.plot(xx, -yy3, 'k--', label=r'$\left|u^{\prime} w^{\prime}\right| = 1  \overline{\left|u^{\prime} w^{\prime}\right|}$',linewidth = 2)
plt.axhline(y=0, xmin=0, xmax=1, c = "k")
plt.axvline(x=0, ymin=0, ymax=1, c = "k")
plt.grid(linestyle = ':')
plt.setp(ax, xlim = (-15, 15), ylim = (-6, 6))
ax.set_xlabel(r"$u'~ [\mathrm{m s^{-1}}$]", fontsize = 20)
ax.set_ylabel(r"$w'~ [\mathrm{m s^{-1}}$]", fontsize = 20)
legend = plt.legend(loc=(0, 0.14), frameon=False,fontsize = 18)
ax.text(y = 4.5,    x = 8, s = "Q1 outward \n interaction", fontweight='heavy', fontsize = 16, c = "k")
ax.text(y = -5.5, x = 8, s = "Q4  sweep", fontweight='heavy', fontsize = 16, c = "k")
ax.text(y = 5,    x = -14, s = "Q2  ejection", fontweight='heavy', fontsize = 16, c = "k")
ax.text(y = -5.5, x = -14, s = "Q3 inward \n interaction", fontweight='heavy', fontsize = 16, c = "k")
plt.tight_layout()
plt.savefig('example.png', dpi = 200)
#%% quadrant plot of W12 OP
fig1 = plt.figure(figsize=(6, 6))
ax = fig1.add_subplot(1,1,1)
h0 = ax.plot( HP , Q_W12_RE2[:, 0], 'k-+', label='Q1',linewidth = 2, ms = 4)
h1 = ax.plot(-HP,  Q_W12_RE2[:, 1], 'k-+', label='Q2',linewidth = 2, ms = 4)
h2 = ax.plot(-HP, -Q_W12_RE2[:, 2], 'k-+', label='Q3',linewidth = 2, ms = 4)
h3 = ax.plot( HP ,-Q_W12_RE2[:, 3], 'k-+', label='Q4',linewidth = 2, ms = 4)

h0 = ax.plot( HP , Q_W12_OP[:, 0], 'r-+', label='Q1',linewidth = 2, ms = 4)
h1 = ax.plot(-HP,  Q_W12_OP[:, 1], 'r-+', label='Q2',linewidth = 2, ms = 4)
h2 = ax.plot(-HP, -Q_W12_OP[:, 2], 'r-+', label='Q3',linewidth = 2, ms = 4)
h3 = ax.plot( HP ,-Q_W12_OP[:, 3], 'r-+', label='Q4',linewidth = 2, ms = 4)
ax.set_title("W12 HF duration fraction", fontsize = 18)

# legend = plt.legend(loc='best', frameon=False,fontsize = 18)
plt.grid(linestyle = ':')
plt.tight_layout()

#%%
HP = np.arange(0, 20, 0.1)
Q_W31_OP = []
for h in HP:
    Q = QH_analysis(W31_fluc_OP, h)
    Q_W31_OP = np.append(Q_W31_OP, Q)
Q_W31_OP = Q_W31_OP.reshape((HP.shape[0], 8))

Q_W31_RE2 = []
for h in HP:
    Q = QH_analysis(W31_fluc_RE2, h)
    Q_W31_RE2 = np.append(Q_W31_RE2, Q)
Q_W31_RE2 = Q_W31_RE2.reshape((HP.shape[0], 8))

#%%

fig1 = plt.figure(figsize=(6, 6))
ax = fig1.add_subplot(1,1,1)
h0 = ax.plot( HP,  Q_W31_RE2[:, 4], 'k-^', label='Q1',linewidth = 2, ms = 4)
h1 = ax.plot(-HP, -Q_W31_RE2[:, 5], 'k-^', label='Q2',linewidth = 2, ms = 4)
h2 = ax.plot(-HP, -Q_W31_RE2[:, 6], 'k-^', label='Q3',linewidth = 2, ms = 4)
h3 = ax.plot( HP,  Q_W31_RE2[:, 7], 'k-^', label='Q4',linewidth = 2, ms = 4)

h0 = ax.plot( HP,  Q_W31_OP[:, 4], 'r-^', label='Q1',linewidth = 2, ms = 4)
h1 = ax.plot(-HP, -Q_W31_OP[:, 5], 'r-^', label='Q2',linewidth = 2, ms = 4)
h2 = ax.plot(-HP, -Q_W31_OP[:, 6], 'r-^', label='Q3',linewidth = 2, ms = 4)
h3 = ax.plot( HP,  Q_W31_OP[:, 7], 'r-^', label='Q4',linewidth = 2, ms = 4)
# plt.yscale('symlog')
ax.set_title("W31 heat flux fraction", fontsize = 18)
# legend = plt.legend(loc='best', frameon=False,fontsize = 18)
plt.grid(linestyle = ':')
plt.tight_layout()
# plt.savefig('output_plots/...')

#%%
fig1 = plt.figure(figsize=(6, 6))
ax = fig1.add_subplot(1,1,1)

h0 = ax.plot( HP , Q_W31_RE2[:, 0], 'k-^', label='Q1',linewidth = 2, ms = 4)
h1 = ax.plot(-HP,  Q_W31_RE2[:, 1], 'k-^', label='Q2',linewidth = 2, ms = 4)
h2 = ax.plot(-HP, -Q_W31_RE2[:, 2], 'k-^', label='Q3',linewidth = 2, ms = 4)
h3 = ax.plot( HP ,-Q_W31_RE2[:, 3], 'k-^', label='Q4',linewidth = 2, ms = 4)

h0 = ax.plot( HP , Q_W31_OP[:, 0], 'r-^', label='Q1',linewidth = 2, ms = 4)
h1 = ax.plot(-HP,  Q_W31_OP[:, 1], 'r-^', label='Q2',linewidth = 2, ms = 4)
h2 = ax.plot(-HP, -Q_W31_OP[:, 2], 'r-^', label='Q3',linewidth = 2, ms = 4)
h3 = ax.plot( HP ,-Q_W31_OP[:, 3], 'r-^', label='Q4',linewidth = 2, ms = 4)

ax.set_title("W31 HF duration fraction", fontsize = 18)
# legend = plt.legend(loc='best', frameon=False,fontsize = 18)
plt.grid(linestyle = ':')
plt.tight_layout()
# plt.savefig('output_plots/...')
#%% report the value tabel

def QH_analysis_uw(W12_fluc_OP, H):
    # thres_H = H * np.std(W12_fluc_OP.w) * np.std(W12_fluc_OP.u)
    thres_H = H * np.mean(np.abs(W12_fluc_OP['uw']))
    # thres_H = H * np.abs(np.mean(W12_fluc_OP['uw']))
    # w>0, u>0
    Q1_W12 = W12_fluc_OP[(W12_fluc_OP["w"]>0)&(W12_fluc_OP["u"]>0)&(np.abs(W12_fluc_OP['uw'])>=thres_H)]
    # w<0, u>0
    Q4_W12 = W12_fluc_OP[(W12_fluc_OP["w"]<0)&(W12_fluc_OP["u"]>0)&(np.abs(W12_fluc_OP['uw'])>=thres_H)]
    # w<0, u<0
    Q3_W12 = W12_fluc_OP[(W12_fluc_OP["w"]<0)&(W12_fluc_OP["u"]<0)&(np.abs(W12_fluc_OP['uw'])>=thres_H)]
    # w>0, u<0
    Q2_W12 = W12_fluc_OP[(W12_fluc_OP["w"]>0)&(W12_fluc_OP["u"]<0)&(np.abs(W12_fluc_OP['uw'])>=thres_H)]

    Q1_tiH = Q1_W12.shape[0] / W12_fluc_OP.shape[0]
    Q2_tiH = Q2_W12.shape[0] / W12_fluc_OP.shape[0]
    Q3_tiH = Q3_W12.shape[0] / W12_fluc_OP.shape[0]
    Q4_tiH = Q4_W12.shape[0] / W12_fluc_OP.shape[0]

    # Q1_SiH = 1 / W12_fluc_OP.shape[0] * Q1_W12["uw"].sum() / (np.std(W12_fluc_OP.w) * np.std(W12_fluc_OP.u))
    # Q2_SiH = 1 / W12_fluc_OP.shape[0] * Q2_W12["uw"].sum() / (np.std(W12_fluc_OP.w) * np.std(W12_fluc_OP.u))
    # Q3_SiH = 1 / W12_fluc_OP.shape[0] * Q3_W12["uw"].sum() / (np.std(W12_fluc_OP.w) * np.std(W12_fluc_OP.u))
    # Q4_SiH = 1 / W12_fluc_OP.shape[0] * Q4_W12["uw"].sum() / (np.std(W12_fluc_OP.w) * np.std(W12_fluc_OP.u))
    
    Q1_SiH = Q1_W12["uw"].sum() / np.abs(W12_fluc_OP['uw']).sum()
    Q2_SiH = Q2_W12["uw"].sum() / np.abs(W12_fluc_OP['uw']).sum()
    Q3_SiH = Q3_W12["uw"].sum() / np.abs(W12_fluc_OP['uw']).sum()
    Q4_SiH = Q4_W12["uw"].sum() / np.abs(W12_fluc_OP['uw']).sum()
    
    return([Q1_tiH, Q2_tiH, Q3_tiH, Q4_tiH, Q1_SiH, Q2_SiH, Q3_SiH, Q4_SiH])
#%%
HP = np.arange(0, 20, 0.1)
Q_W12_OP_uw = []
for h in HP:
    Q = QH_analysis_uw(W12_fluc_OP, h)
    Q_W12_OP_uw = np.append(Q_W12_OP_uw, Q)
Q_W12_OP_uw = Q_W12_OP_uw.reshape((HP.shape[0], 8))

Q_W12_RE2_uw = []
for h in HP:
    Q = QH_analysis_uw(W12_fluc_RE2, h)
    Q_W12_RE2_uw = np.append(Q_W12_RE2_uw, Q)
Q_W12_RE2_uw = Q_W12_RE2_uw.reshape((HP.shape[0], 8))


#%% quadrant plot of W12 OP
fig1 = plt.figure(figsize=(6, 6))
ax = fig1.add_subplot(1,1,1)

h0 = ax.plot( HP,  Q_W12_RE2_uw[:, 4], 'k-+', label='Q1',linewidth = 2, ms = 4)
h1 = ax.plot(-HP, -Q_W12_RE2_uw[:, 5], 'k-+', label='Q2',linewidth = 2, ms = 4)
h2 = ax.plot(-HP, -Q_W12_RE2_uw[:, 6], 'k-+', label='Q3',linewidth = 2, ms = 4)
h3 = ax.plot( HP,  Q_W12_RE2_uw[:, 7], 'k-+', label='Q4',linewidth = 2, ms = 4)

h0 = ax.plot( HP,  Q_W12_OP_uw[:, 4], 'r-+', label='Q1',linewidth = 2, ms = 4)
h1 = ax.plot(-HP, -Q_W12_OP_uw[:, 5], 'r-+', label='Q2',linewidth = 2, ms = 4)
h2 = ax.plot(-HP, -Q_W12_OP_uw[:, 6], 'r-+', label='Q3',linewidth = 2, ms = 4)
h3 = ax.plot( HP,  Q_W12_OP_uw[:, 7], 'r-+', label='Q4',linewidth = 2, ms = 4)
ax.set_title("W12 stress fraction", fontsize = 18)
# plt.yscale('symlog')
# plt.xscale('symlog')
# legend = plt.legend(loc='best', frameon=False,fontsize = 18)
plt.grid(linestyle = ':')
plt.tight_layout()
# plt.savefig('output_plots/...')

#%% quadrant plot of W12 OP
fig1 = plt.figure(figsize=(6, 6))
ax = fig1.add_subplot(1,1,1)
h0 = ax.plot( HP , Q_W12_RE2_uw[:, 0], 'k-+', label='Q1',linewidth = 2, ms = 4)
h1 = ax.plot(-HP,  Q_W12_RE2_uw[:, 1], 'k-+', label='Q2',linewidth = 2, ms = 4)
h2 = ax.plot(-HP, -Q_W12_RE2_uw[:, 2], 'k-+', label='Q3',linewidth = 2, ms = 4)
h3 = ax.plot( HP ,-Q_W12_RE2_uw[:, 3], 'k-+', label='Q4',linewidth = 2, ms = 4)

h0 = ax.plot( HP , Q_W12_OP_uw[:, 0], 'r-+', label='Q1',linewidth = 2, ms = 4)
h1 = ax.plot(-HP,  Q_W12_OP_uw[:, 1], 'r-+', label='Q2',linewidth = 2, ms = 4)
h2 = ax.plot(-HP, -Q_W12_OP_uw[:, 2], 'r-+', label='Q3',linewidth = 2, ms = 4)
h3 = ax.plot( HP ,-Q_W12_OP_uw[:, 3], 'r-+', label='Q4',linewidth = 2, ms = 4)
ax.set_title("W12 stress duration fraction", fontsize = 18)

# legend = plt.legend(loc='best', frameon=False,fontsize = 18)
plt.grid(linestyle = ':')
plt.tight_layout()

#%%
HP = np.arange(0, 20, 0.1)
Q_W31_OP_uw = []
for h in HP:
    Q = QH_analysis_uw(W31_fluc_OP, h)
    Q_W31_OP_uw = np.append(Q_W31_OP_uw, Q)
Q_W31_OP_uw = Q_W31_OP_uw.reshape((HP.shape[0], 8))


Q_W31_RE2_uw = []
for h in HP:
    Q = QH_analysis_uw(W31_fluc_RE2, h)
    Q_W31_RE2_uw = np.append(Q_W31_RE2_uw, Q)
Q_W31_RE2_uw = Q_W31_RE2_uw.reshape((HP.shape[0], 8))
#%%

fig1 = plt.figure(figsize=(6, 6))
ax = fig1.add_subplot(1,1,1)
h0 = ax.plot( HP,  Q_W31_RE2_uw[:, 4], 'k-^', label='Q1',linewidth = 2, ms = 4)
h1 = ax.plot(-HP, -Q_W31_RE2_uw[:, 5], 'k-^', label='Q2',linewidth = 2, ms = 4)
h2 = ax.plot(-HP, -Q_W31_RE2_uw[:, 6], 'k-^', label='Q3',linewidth = 2, ms = 4)
h3 = ax.plot( HP,  Q_W31_RE2_uw[:, 7], 'k-^', label='Q4',linewidth = 2, ms = 4)

h0 = ax.plot( HP,  Q_W31_OP_uw[:, 4], 'r-^', label='Q1',linewidth = 2, ms = 4)
h1 = ax.plot(-HP, -Q_W31_OP_uw[:, 5], 'r-^', label='Q2',linewidth = 2, ms = 4)
h2 = ax.plot(-HP, -Q_W31_OP_uw[:, 6], 'r-^', label='Q3',linewidth = 2, ms = 4)
h3 = ax.plot( HP,  Q_W31_OP_uw[:, 7], 'r-^', label='Q4',linewidth = 2, ms = 4)
# plt.yscale('symlog')
# plt.xscale('symlog')
ax.set_title("W31 stress fraction", fontsize = 18)
# legend = plt.legend(loc='best', frameon=False,fontsize = 18)
plt.grid(linestyle = ':')
plt.tight_layout()
# plt.savefig('output_plots/...')

#%%
fig1 = plt.figure(figsize=(6, 6))
ax = fig1.add_subplot(1,1,1)

h0 = ax.plot( HP , Q_W31_RE2_uw[:, 0], 'k-^', label='Q1',linewidth = 2, ms = 4)
h1 = ax.plot(-HP,  Q_W31_RE2_uw[:, 1], 'k-^', label='Q2',linewidth = 2, ms = 4)
h2 = ax.plot(-HP, -Q_W31_RE2_uw[:, 2], 'k-^', label='Q3',linewidth = 2, ms = 4)
h3 = ax.plot( HP ,-Q_W31_RE2_uw[:, 3], 'k-^', label='Q4',linewidth = 2, ms = 4)

h0 = ax.plot( HP , Q_W31_OP_uw[:, 0], 'r-^', label='Q1',linewidth = 2, ms = 4)
h1 = ax.plot(-HP,  Q_W31_OP_uw[:, 1], 'r-^', label='Q2',linewidth = 2, ms = 4)
h2 = ax.plot(-HP, -Q_W31_OP_uw[:, 2], 'r-^', label='Q3',linewidth = 2, ms = 4)
h3 = ax.plot( HP ,-Q_W31_OP_uw[:, 3], 'r-^', label='Q4',linewidth = 2, ms = 4)

ax.set_title("W31 stress duration fraction", fontsize = 18)
# legend = plt.legend(loc='best', frameon=False,fontsize = 18)
plt.grid(linestyle = ':')
plt.tight_layout()
# plt.savefig('output_plots/...')

#%%
# --------------------------------------------------------- #
# get collection of four plots   sharex = True
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize = (10, 10), sharey = True)
# ax = fig1.add_subplot(1,1,1)
h0 = ax1.plot( HP,  Q_W12_RE2[:, 4], 'k-+', linewidth = 2, ms = 4)
h1 = ax1.plot(-HP, -Q_W12_RE2[:, 5], 'k-+', linewidth = 2, ms = 4)
h2 = ax1.plot(-HP, -Q_W12_RE2[:, 6], 'k-+', linewidth = 2, ms = 4)
h3 = ax1.plot( HP,  Q_W12_RE2[:, 7], 'k-+', linewidth = 2, ms = 4)
h0 = ax1.plot( HP,  Q_W12_OP[:, 4], 'r-+', linewidth = 2, ms = 4)
h1 = ax1.plot(-HP, -Q_W12_OP[:, 5], 'r-+', linewidth = 2, ms = 4)
h2 = ax1.plot(-HP, -Q_W12_OP[:, 6], 'r-+', linewidth = 2, ms = 4)
h3 = ax1.plot( HP,  Q_W12_OP[:, 7], 'r-+', linewidth = 2, ms = 4)

h0 = ax2.plot( HP,  Q_W31_RE2[:, 4], 'k-^', label='off', linewidth = 2, ms = 4)
h1 = ax2.plot(-HP, -Q_W31_RE2[:, 5], 'k-^',linewidth = 2, ms = 4)
h2 = ax2.plot(-HP, -Q_W31_RE2[:, 6], 'k-^',linewidth = 2, ms = 4)
h3 = ax2.plot( HP,  Q_W31_RE2[:, 7], 'k-^',linewidth = 2, ms = 4)
h0 = ax2.plot( HP,  Q_W31_OP[:, 4], 'r-^', label='on', linewidth = 2, ms = 4)
h1 = ax2.plot(-HP, -Q_W31_OP[:, 5], 'r-^',linewidth = 2, ms = 4)
h2 = ax2.plot(-HP, -Q_W31_OP[:, 6], 'r-^',linewidth = 2, ms = 4)
h3 = ax2.plot( HP,  Q_W31_OP[:, 7], 'r-^',linewidth = 2, ms = 4)

h0 = ax3.plot( HP,  Q_W12_RE2_uw[:, 4], 'k-+',linewidth = 2, ms = 4)
h1 = ax3.plot(-HP, -Q_W12_RE2_uw[:, 5], 'k-+',linewidth = 2, ms = 4)
h2 = ax3.plot(-HP, -Q_W12_RE2_uw[:, 6], 'k-+',linewidth = 2, ms = 4)
h3 = ax3.plot( HP,  Q_W12_RE2_uw[:, 7], 'k-+',linewidth = 2, ms = 4)
h0 = ax3.plot( HP,  Q_W12_OP_uw[:, 4], 'r-+',linewidth = 2, ms = 4)
h1 = ax3.plot(-HP, -Q_W12_OP_uw[:, 5], 'r-+',linewidth = 2, ms = 4)
h2 = ax3.plot(-HP, -Q_W12_OP_uw[:, 6], 'r-+',linewidth = 2, ms = 4)
h3 = ax3.plot( HP,  Q_W12_OP_uw[:, 7], 'r-+',linewidth = 2, ms = 4)

h0 = ax4.plot( HP,  Q_W31_RE2_uw[:, 4], 'k-^',linewidth = 2, ms = 4)
h1 = ax4.plot(-HP, -Q_W31_RE2_uw[:, 5], 'k-^',linewidth = 2, ms = 4)
h2 = ax4.plot(-HP, -Q_W31_RE2_uw[:, 6], 'k-^',linewidth = 2, ms = 4)
h3 = ax4.plot( HP,  Q_W31_RE2_uw[:, 7], 'k-^',linewidth = 2, ms = 4)
h0 = ax4.plot( HP,  Q_W31_OP_uw[:, 4], 'r-^', linewidth = 2, ms = 4)
h1 = ax4.plot(-HP, -Q_W31_OP_uw[:, 5], 'r-^', linewidth = 2, ms = 4)
h2 = ax4.plot(-HP, -Q_W31_OP_uw[:, 6], 'r-^', linewidth = 2, ms = 4)
h3 = ax4.plot( HP,  Q_W31_OP_uw[:, 7], 'r-^', linewidth = 2, ms = 4)
# ax.set_title("W31 stress duration fraction", fontsize = 18)
# legend = plt.legend(loc='best', frameon=False,fontsize = 18)
ax1.grid(linestyle = ':')
ax2.grid(linestyle = ':')
ax3.grid(linestyle = ':')
ax4.grid(linestyle = ':')
ax1.text(x = 14, y = -0.66, s = r"\textbf{a)A1}", fontweight='heavy', fontsize = 20, c = "k")
ax2.text(x = 14, y = -0.66, s = r"\textbf{b)A2}", fontweight='heavy', fontsize = 20, c = "k")
ax3.text(x = 14, y = -0.66, s = r"\textbf{c)A1}", fontweight='heavy', fontsize = 20, c = "k")
ax4.text(x = 14, y = -0.66, s = r"\textbf{d)A2}", fontweight='heavy', fontsize = 20, c = "k")

ax1.text(x = 6, y = 0.20, s = "Q1 warm \n updraft", fontweight='heavy', fontsize = 20, c = "k")
ax1.text(x =-15, y = 0.20, s = "Q2 cool \n updraft", fontweight='heavy', fontsize = 20, c = "k")
ax1.text(x = 8, y =-0.45, s = "Q4 warm \n downdraft", fontweight='heavy', fontsize = 20, c = "k")
ax1.text(x =-15, y =-0.45, s = "Q3 cool \n downdraft", fontweight='heavy', fontsize = 20, c = "k")

ax3.text(x = 6, y = 0.20, s = "Q1 outward \n interaction", fontweight='heavy', fontsize = 20, c = "k")
ax3.text(x =-15, y = 0.20, s = "Q2 \n ejection", fontweight='heavy', fontsize = 20, c = "k")
ax3.text(x = 8, y =-0.45, s = "Q4 \n sweep", fontweight='heavy', fontsize = 20, c = "k")
ax3.text(x =-15, y =-0.45, s = "Q3 inward \n interaction", fontweight='heavy', fontsize = 20, c = "k")

xxlabels = [item.get_text() for item in ax3.get_xticklabels()]
xxlabels = ['$\\mathdefault{20}$',
 '$\\mathdefault{10}$',
 '$\\mathdefault{0}$',
 '$\\mathdefault{10}$',
 '$\\mathdefault{20}$']
ax3.set_xticklabels(xxlabels)
yylabels = [item.get_text() for item in ax3.get_xticklabels()]
yylabels = ['$\\mathdefault{0.8}$',
 '$\\mathdefault{0.6}$',
 '$\\mathdefault{0.4}$',
 '$\\mathdefault{0.2}$',
 '$\\mathdefault{0.0}$',
 '$\\mathdefault{0.2}$',
 '$\\mathdefault{0.4}$']
ax3.set_yticklabels(yylabels)
plt.setp((ax1, ax2), ylim = (-0.7, 0.4), xlim = (-20, 20))
plt.setp((ax3, ax4), ylim = (-0.7, 0.4), xlim = (-20, 20))
ax1.set_ylabel(r"Heat flux \textit{value fraction}",fontsize = 20)
ax3.set_ylabel(r"Momentum flux \textit{value fraction}",fontsize = 20)
ax3.set_xlabel("Hole size H",fontsize = 20)
ax4.set_xlabel("Hole size H",fontsize = 20)
# legend = plt.legend(loc='best', frameon=False,fontsize = 20)
ax2.legend(loc='upper left', frameon=False,fontsize = 20)
plt.tight_layout()
plt.savefig("QH_flux.pdf")

# %% ---------------------------------------------------------- #
#   -- report quadrant hole anslysis key parameter --- # 
# % ----------------------------------------------------------- #
# correlation between two variable
wT_A1 = [np.corrcoef(W12_fluc_OP.w,  W12_fluc_OP.Ts)[0,1],
         np.corrcoef(W12_fluc_RE2.w, W12_fluc_RE2.Ts)[0,1],
         np.corrcoef(W31_fluc_OP.w,  W31_fluc_OP.Ts)[0,1],
         np.corrcoef(W31_fluc_RE2.w, W31_fluc_RE2.Ts)[0,1]]

# find where half momentum transfer occurs 
def find_nearest_idx(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

# H' is the hole size above which half of flux occurs
wT_A20 = [find_nearest_idx(np.sum(np.abs(Q_W12_OP[:,4:]), axis = 1), 0.5),
      find_nearest_idx(np.sum(np.abs(Q_W12_RE2[:,4:]), axis = 1), 0.5),
      find_nearest_idx(np.sum(np.abs(Q_W31_OP[:,4:]), axis = 1), 0.5),
      find_nearest_idx(np.sum(np.abs(Q_W31_RE2[:,4:]), axis = 1), 0.5)]
wT_A2 = np.array(wT_A20) * 0.1

# sum(tiH')
wT_A3 = [np.sum(Q_W12_OP[:,:4], axis = 1)[wT_A20[0]],
      np.sum(Q_W12_RE2[:,:4], axis = 1)[wT_A20[1]],
      np.sum(Q_W31_OP[:,:4], axis = 1)[wT_A20[2]],
      np.sum(Q_W31_RE2[:,:4], axis = 1)[wT_A20[3]]]

# F30 / F10
# wT_A4 = [Q_W12_OP[0, 6]/Q_W12_OP[0, 4],
#          Q_W12_RE2[0, 6]/Q_W12_RE2[0, 4],
#          Q_W31_OP[0, 6]/Q_W31_OP[0, 4],
#          Q_W31_RE2[0, 6]/Q_W31_RE2[0, 4]]
# F40 / F20
wT_A4 = [ Q_W12_OP[0, 7]/ Q_W12_OP[0, 5],
         Q_W12_RE2[0, 7]/Q_W12_RE2[0, 5],
          Q_W31_OP[0, 7]/ Q_W31_OP[0, 5],
         Q_W31_RE2[0, 7]/Q_W31_RE2[0, 5]]

# (F20 + F40) / (F30 + F10)
wT_A5 = [(Q_W12_OP[0, 5] + Q_W12_OP[0, 7]) / (Q_W12_OP[0, 6] + Q_W12_OP[0, 4]),
         (Q_W12_RE2[0, 5] + Q_W12_RE2[0, 7]) / (Q_W12_RE2[0, 6] +Q_W12_RE2[0, 4]),
         (Q_W31_OP[0, 5] + Q_W31_OP[0, 7]) / (Q_W31_OP[0, 6] + Q_W31_OP[0, 4]),
         (Q_W31_RE2[0, 5] + Q_W31_RE2[0, 7]) / (Q_W31_RE2[0, 6] +Q_W31_RE2[0, 4])]

# F3h' / F1h'
# wT_A6 = [Q_W12_OP[wT_A20[0], 6]/Q_W12_OP[wT_A20[0], 4],
#       Q_W12_RE2[wT_A20[1], 6]/Q_W12_RE2[wT_A20[1], 4],
#       Q_W31_OP[wT_A20[2], 6]/Q_W31_OP[wT_A20[2], 4],
#       Q_W31_RE2[wT_A20[3], 6]/Q_W31_RE2[wT_A20[3], 4]]

# F4h' / F2h'
wT_A6 = [Q_W12_OP[wT_A20[0], 7]/ Q_W12_OP[wT_A20[0], 5],
        Q_W12_RE2[wT_A20[1], 7]/Q_W12_RE2[wT_A20[1], 5],
         Q_W31_OP[wT_A20[2], 7]/ Q_W31_OP[wT_A20[2], 5],
        Q_W31_RE2[wT_A20[3], 7]/Q_W31_RE2[wT_A20[3], 5]]

# (F2h + F4h) / (F3h + F1h)
wT_A7 = [(Q_W12_OP[wT_A20[0], 5] + Q_W12_OP[wT_A20[0], 7]) / (Q_W12_OP[wT_A20[0], 6] + Q_W12_OP[wT_A20[0], 4]),
         (Q_W12_RE2[wT_A20[1], 5] + Q_W12_RE2[wT_A20[1], 7]) / (Q_W12_RE2[wT_A20[1], 6] +Q_W12_RE2[wT_A20[1], 4]),
         (Q_W31_OP[wT_A20[2], 5] + Q_W31_OP[wT_A20[2], 7]) / (Q_W31_OP[wT_A20[2], 6] + Q_W31_OP[wT_A20[2], 4]),
         (Q_W31_RE2[wT_A20[3], 5] + Q_W31_RE2[wT_A20[3], 7]) / (Q_W31_RE2[wT_A20[3], 6] +Q_W31_RE2[wT_A20[3], 4])]
#%%
# --------------------------------------------------------- #
uw_A1 = [np.corrcoef(W12_fluc_OP.w,  W12_fluc_OP.u)[0,1],
      np.corrcoef(W12_fluc_RE2.w, W12_fluc_RE2.u)[0,1],
      np.corrcoef(W31_fluc_OP.w,  W31_fluc_OP.u)[0,1],
      np.corrcoef(W31_fluc_RE2.w, W31_fluc_RE2.u)[0,1]]

# H' is the hole size above which half of flux occurs
uw_A20 = [find_nearest_idx(np.sum(np.abs(Q_W12_OP_uw[:,4:]), axis = 1), 0.5),
      find_nearest_idx(np.sum(np.abs(Q_W12_RE2_uw[:,4:]), axis = 1), 0.5),
      find_nearest_idx(np.sum(np.abs(Q_W31_OP_uw[:,4:]), axis = 1), 0.5),
      find_nearest_idx(np.sum(np.abs(Q_W31_RE2_uw[:,4:]), axis = 1), 0.5)]
uw_A2 = np.array(uw_A20) * 0.1

# sum(tiH')
uw_A3 = [np.sum(Q_W12_OP_uw[:,:4], axis = 1)[uw_A20[0]],
      np.sum(Q_W12_RE2_uw[:,:4], axis = 1)[uw_A20[1]],
      np.sum(Q_W31_OP_uw[:,:4], axis = 1)[uw_A20[2]],
      np.sum(Q_W31_RE2_uw[:,:4], axis = 1)[uw_A20[3]]]

# F30 / F10
# uw_A4 = [Q_W12_OP_uw[0, 6]/Q_W12_OP_uw[0, 4],
#          Q_W12_RE2_uw[0, 6]/Q_W12_RE2_uw[0, 4],
#          Q_W31_OP_uw[0, 6]/Q_W31_OP_uw[0, 4],
#          Q_W31_RE2_uw[0, 6]/Q_W31_RE2_uw[0, 4]]
# F40 / F20
uw_A4 = [Q_W12_OP_uw[0, 7]/ Q_W12_OP_uw[0, 5],
        Q_W12_RE2_uw[0, 7]/Q_W12_RE2_uw[0, 5],
         Q_W31_OP_uw[0, 7]/ Q_W31_OP_uw[0, 5],
        Q_W31_RE2_uw[0, 7]/Q_W31_RE2_uw[0, 5]]


# (F20 + F40) / (F30 + F10)
uw_A5 = [(Q_W12_OP_uw[0, 5] + Q_W12_OP_uw[0, 7]) / (Q_W12_OP_uw[0, 6] + Q_W12_OP_uw[0, 4]),
         (Q_W12_RE2_uw[0, 5] + Q_W12_RE2_uw[0, 7]) / (Q_W12_RE2_uw[0, 6] +Q_W12_RE2_uw[0, 4]),
         (Q_W31_OP_uw[0, 5] + Q_W31_OP_uw[0, 7]) / (Q_W31_OP_uw[0, 6] + Q_W31_OP_uw[0, 4]),
         (Q_W31_RE2_uw[0, 5] + Q_W31_RE2_uw[0, 7]) / (Q_W31_RE2_uw[0, 6] +Q_W31_RE2_uw[0, 4])]

# F3h' / F1h'
# uw_A6 = [Q_W12_OP_uw[uw_A20[0], 6]/Q_W12_OP_uw[uw_A20[0], 4],
#       Q_W12_RE2_uw[uw_A20[1], 6]/Q_W12_RE2_uw[uw_A20[1], 4],
#       Q_W31_OP_uw[uw_A20[2], 6]/Q_W31_OP_uw[uw_A20[2], 4],
#       Q_W31_RE2_uw[uw_A20[3], 6]/Q_W31_RE2_uw[uw_A20[3], 4]]
# F4h' / F0h'
uw_A6 = [Q_W12_OP_uw[uw_A20[0], 7]/ Q_W12_OP_uw[uw_A20[0], 5],
        Q_W12_RE2_uw[uw_A20[1], 7]/Q_W12_RE2_uw[uw_A20[1], 5],
         Q_W31_OP_uw[uw_A20[2], 7]/ Q_W31_OP_uw[uw_A20[2], 5],
        Q_W31_RE2_uw[uw_A20[3], 7]/Q_W31_RE2_uw[uw_A20[3], 5]]
# (F2h + F4h) / (F3h + F1h)
uw_A7 = [(Q_W12_OP_uw[uw_A20[0], 5] + Q_W12_OP_uw[uw_A20[0], 7]) / (Q_W12_OP_uw[uw_A20[0], 6] + Q_W12_OP_uw[uw_A20[0], 4]),
         (Q_W12_RE2_uw[uw_A20[1], 5] + Q_W12_RE2_uw[uw_A20[1], 7]) / (Q_W12_RE2_uw[uw_A20[1], 6] +Q_W12_RE2_uw[uw_A20[1], 4]),
         (Q_W31_OP_uw[uw_A20[2], 5] + Q_W31_OP_uw[uw_A20[2], 7]) / (Q_W31_OP_uw[uw_A20[2], 6] + Q_W31_OP_uw[uw_A20[2], 4]),
         (Q_W31_RE2_uw[uw_A20[3], 5] + Q_W31_RE2_uw[uw_A20[3], 7]) / (Q_W31_RE2_uw[uw_A20[3], 6] +Q_W31_RE2_uw[uw_A20[3], 4])]

#%% 











































# #%% W12 quadrant analysis of the flux during OP

# # ----- analysis on w'T' ----- 
# fig1 = plt.figure(figsize=(6, 6))
# ax = fig1.add_subplot(1,1,1)
# ax.scatter(W12_fluc_OP.w, W12_fluc_OP.Ts, s = np.abs(W12_WT_OP)*3,  marker = ".")
# ax.set_ylabel("T' [$^\circ C$]", fontsize = 18)
# ax.set_xlabel("w' [$m s^{-1}$]", fontsize = 18)
# ax.set_title("W12 OP", fontsize = 18)
# plt.grid(linestyle = ':')
# plt.tight_layout()

# #%% W12 quadrant analysis of the flux during RE2
# fig1 = plt.figure(figsize=(6, 6))
# ax = fig1.add_subplot(1,1,1)
# plt.scatter(W12_fluc_RE2.w, W12_fluc_RE2.Ts, s = np.abs(W12_WT_RE2)*30,  marker = ".")
# # plt.scatter(W12_fluc_RE2.w, W12_fluc_RE2.Ts, s = 1,  marker = ".")
# ax.set_ylabel("T' [$^\circ C$]", fontsize = 18)
# ax.set_xlabel("w' [$m s^{-1}$]", fontsize = 18)
# ax.set_title("W12 RE2", fontsize = 18)
# plt.grid(linestyle = ':')
# plt.tight_layout()

# #%% W31 quadrant analysis of the flux during OP
# fig1 = plt.figure(figsize=(6, 6))
# ax = fig1.add_subplot(1,1,1)
# plt.scatter(W31_fluc_OP.w, W31_fluc_OP.Ts, s = 1,  marker = ".")
# ax.set_ylabel("T' [$^\circ C$]", fontsize = 18)
# ax.set_xlabel("w' [$m s^{-1}$]", fontsize = 18)
# ax.set_title("W31 OP", fontsize = 18)
# plt.grid(linestyle = ':')
# plt.tight_layout()

# #%% W31 quadrant analysis of the flux during RE2
# fig1 = plt.figure(figsize=(6, 6))
# ax = fig1.add_subplot(1,1,1)
# plt.scatter(W31_fluc_RE2.w, W31_fluc_RE2.Ts, s = 1,  marker = ".")
# ax.set_ylabel("T' [$^\circ C$]", fontsize = 18)
# ax.set_xlabel("w' [$m s^{-1}$]", fontsize = 18)
# ax.set_title("W31 RE2", fontsize = 18)
# plt.grid(linestyle = ':')
# plt.tight_layout()

# #%%%
# # ----- analysis on u'w' -----
# #%% W12 quadrant analysis of the flux during OP
# fig1 = plt.figure(figsize=(6, 6))
# ax = fig1.add_subplot(1,1,1)
# ax.scatter(W12_fluc_OP.w, W12_fluc_OP.u, s = 1,  marker = ".")
# ax.set_ylabel("u' [$^\circ C$]", fontsize = 18)
# ax.set_xlabel("w' [$m s^{-1}$]", fontsize = 18)
# ax.set_title("W12 OP ", fontsize = 18)
# plt.grid(linestyle = ':')
# plt.tight_layout()

# #%% W12 quadrant analysis of the flux during RE2
# fig1 = plt.figure(figsize=(6, 6))
# ax = fig1.add_subplot(1,1,1)
# plt.scatter(W12_fluc_RE2.w, W12_fluc_RE2.u, s = 1,  marker = ".")
# ax.set_ylabel("u' [$^\circ C$]", fontsize = 18)
# ax.set_xlabel("w' [$m s^{-1}$]", fontsize = 18)
# ax.set_title("W12 RE2", fontsize = 18)
# plt.grid(linestyle = ':')
# plt.tight_layout()

# #%% W31 quadrant analysis of the flux during OP
# fig1 = plt.figure(figsize=(6, 6))
# ax = fig1.add_subplot(1,1,1)
# plt.scatter(W31_fluc_OP.w, W31_fluc_OP.u, s = 1,  marker = ".")
# ax.set_ylabel("u' [$^\circ C$]", fontsize = 18)
# ax.set_xlabel("w' [$m s^{-1}$]", fontsize = 18)
# ax.set_title("W31 OP", fontsize = 18)
# plt.grid(linestyle = ':')
# plt.tight_layout()

# #%% W31 quadrant analysis of the flux during RE2
# fig1 = plt.figure(figsize=(6, 6))
# ax = fig1.add_subplot(1,1,1)
# plt.scatter(W31_fluc_RE2.w, W31_fluc_RE2.u, s = 1,  marker = ".")
# ax.set_ylabel("u' [$^\circ C$]", fontsize = 18)
# ax.set_xlabel("w' [$m s^{-1}$]", fontsize = 18)
# ax.set_title("W31 RE2", fontsize = 18)
# plt.grid(linestyle = ':')
# plt.tight_layout()




# #%%
# fig1 = plt.figure(figsize=(6, 6))
# ax = fig1.add_subplot(1,1,1)
# h1 = plt.plot(W12_WT_OP.index, W12_WT_OP, 'k--', label='flux',linewidth = 2)

# legend = plt.legend(loc='best', frameon=False,fontsize = 18)
# plt.grid(linestyle = ':')
# plt.tight_layout()
# #%%


# #%%


















