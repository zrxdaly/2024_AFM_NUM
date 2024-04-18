#%% 
# import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import glob
import os
# plt.rc('font', family='serif')
plt.rc('xtick', labelsize='20')
plt.rc('ytick', labelsize='20')
plt.rc('text', usetex=True)

import scienceplots
plt.style.use(['science', 'grid', 'scatter'])
#%%
file_dir = os.path.dirname(os.path.realpath(__file__))
idx = 3
heig = str(np.arange(5, 105, 10)[idx]).zfill(3)
files_BW = sorted(glob.glob(file_dir + "/HS/PAV_UW%s_*.npy"%heig))

#%%
flux_BW_OP = np.load(files_BW[0])
# flux_UW = np.load(files_UW[3])
#%%
fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(1,1,1)
ax.set_aspect(aspect=1)
con = plt.contourf(np.mean(flux_BW_OP, axis = 0), levels = np.arange(-0.8,0.8,0.1), cmap = 'gnuplot2')
# con = plt.contourf(flux_BW[1,:], levels = np.arange(-5,0.5,0.1), cmap = 'gnuplot2')
ax.set_title('', fontsize=18)
cbar = fig.colorbar(con)
cbar.ax.set_ylabel('flux',fontsize = 18)
plt.tight_layout()
# %%
for i in range(flux_BW_OP.shape[0]):
    tt = str(i).zfill(2)
    with plt.ioff():
        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(1,1,1)
        ax.set_aspect(aspect=1)
        con = plt.contourf(flux_BW_OP[i,:], levels = np.arange(-0.5,0.5,0.1), cmap = 'gnuplot2')
        ax.set_title('', fontsize=18)
        cbar = fig.colorbar(con)
        cbar.ax.set_ylabel('flux',fontsize = 18)
        plt.tight_layout()
        ax.set_xlim([1024-500,1024+500])
        ax.set_ylim([0,1000])
        plt.savefig(file_dir + "/HS/UW_plot/%s.png"%tt)
# %%
