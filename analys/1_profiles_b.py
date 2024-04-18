# %%
import numpy as np
import glob
import os 

import matplotlib.pyplot as plt
plt.rc('font', family='serif')
plt.rc('xtick', labelsize='20')
plt.rc('ytick', labelsize='20')
plt.rc('text', usetex=True)

# %% # --------------------------------------------------------- #
file_dir = os.path.dirname(os.path.realpath(__file__))
Bfiles = sorted(glob.glob(file_dir + "/W12/W12_Bt=*x=528"))
BSlice = np.fromfile(Bfiles[0], dtype = float).reshape(629, 11, 11)
for i in range(1, len(Bfiles)-1):
    print(Bfiles[i])
    BSlice0 = np.fromfile(Bfiles[i], dtype = float).reshape(600, 11, 11)
    BSlice = np.append(BSlice, BSlice0, axis = 0)
BSlice0 = np.fromfile(Bfiles[len(Bfiles)-1], dtype = float).reshape(1, 11, 11)
BSlice = np.append(BSlice, BSlice0, axis = 0)
BSlice = BSlice/9.81*273.15
#%% check the temp profile
B_RE0 = np.mean(BSlice[0:100, :, :], axis = 0)
B_RE1 = np.mean(BSlice[100:200, :, :], axis = 0)
B_RE2 = np.mean(BSlice[200:300, :, :], axis = 0)
B_RE3 = np.mean(BSlice[300:400, :, :], axis = 0)
B_OP = np.mean(BSlice[600:1500, :, :], axis = 0)
idx = 4

fig1 = plt.figure(figsize=(6, 9))
ax = fig1.add_subplot(1,1,1)
# for i in range(0, BSlice.shape[0], 5):
#     h1 = plt.plot(BSlice[i, :,5], np.arange(11), '--o', linewidth = 2)
plt.plot(B_OP[:,idx], np.arange(11), 'r-', linewidth = 2)
plt.plot(B_RE0[:,idx], np.arange(11), 'b-', linewidth = 2)
# plt.plot(B_RE1[:,idx], np.arange(11), 'b-', linewidth = 2)
# plt.plot(B_RE2[:,idx], np.arange(11), 'b-', linewidth = 2)
plt.plot(B_RE3[:,idx], np.arange(11), 'b--', linewidth = 2)
plt.xlim([-1,10])
plt.ylim([0,10])
plt.grid(linestyle = ':')
plt.tight_layout()
#%% check temp time series on canopy top
fig1 = plt.figure(figsize=(12, 4))
ax = fig1.add_subplot(1,1,1)
plt.plot(BSlice[:, 1, 3], '--o', linewidth = 2)
# plt.xlim([-1,10])
# plt.ylim([0,10])
plt.grid(linestyle = ':')
plt.tight_layout()
#%% # --------------------------------------------------------- #
# check the velocity profile
Ufiles = sorted(glob.glob(file_dir + "/W12/W12_Ut=*x=528"))
USlice = np.fromfile(Ufiles[0], dtype = float).reshape(629, 11, 11)
for i in range(1, len(Ufiles)-1):
    print(Ufiles[i])
    USlice0 = np.fromfile(Ufiles[i], dtype = float).reshape(600, 11, 11)
    USlice = np.append(USlice, USlice0, axis = 0)
USlice0 = np.fromfile(Ufiles[len(Ufiles)-1], dtype = float).reshape(1, 11, 11)
USlice = np.append(USlice, USlice0, axis = 0)
USlice = USlice
#%%
fig1 = plt.figure(figsize=(6, 9))
ax = fig1.add_subplot(1,1,1)
for i in range(1, USlice.shape[0], 1):
    h1 = plt.plot(USlice[i, :,5], np.arange(11), '--o', linewidth = 2)
# plt.xlim([-1,10])
# plt.ylim([0,10])
plt.grid(linestyle = ':')
plt.tight_layout()
#%% check U time series on canopy top
fig1 = plt.figure(figsize=(9, 4))
ax = fig1.add_subplot(1,1,1)
plt.plot(USlice[:, 2, 4], '--o', linewidth = 2)
# plt.xlim([-1,10])
# plt.ylim([0,10])
plt.grid(linestyle = ':')
plt.tight_layout()
#%% # --------------------------------------------------------- #
# check the Tv and Ta profile
TVfiles = sorted(glob.glob(file_dir + "/W12/W12_TVt=*x=528"))
TVSlice = np.fromfile(TVfiles[0], dtype = float).reshape(629, 11, 11)
for i in range(1, len(TVfiles)-1):
    print(TVfiles[i])
    TVSlice0 = np.fromfile(TVfiles[i], dtype = float).reshape(600, 11, 11)
    TVSlice = np.append(TVSlice, TVSlice0, axis = 0)
TVSlice0 = np.fromfile(TVfiles[len(TVfiles)-1], dtype = float).reshape(1, 11, 11)
TVSlice = np.append(TVSlice, TVSlice0, axis = 0)
TVSlice = TVSlice - 273.15
#%%
fig1 = plt.figure(figsize=(6, 9))
ax = fig1.add_subplot(1,1,1)
for i in range(1, TVSlice.shape[0], 500):
    h1 = plt.plot(TVSlice[i, :,3], np.arange(11), '--o', linewidth = 2)
plt.grid(linestyle = ':')
plt.tight_layout()
#%% check TV and Ta time series on canopy top
idx = 3
fig1 = plt.figure(figsize=(9, 4))
ax = fig1.add_subplot(1,1,1)
plt.plot(TVSlice[:, 2, idx], '--o', linewidth = 1)
plt.plot(BSlice[:, 2, idx], '-', linewidth = 1)
plt.plot(BSlice[:, 2, idx] - TVSlice[:, 2, idx], '+', linewidth = 1)
# plt.xlim([-1,10])
# plt.ylim([0,10])
plt.grid(linestyle = ':')
plt.tight_layout()

# %%
