#%% 
# import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import glob
import os
# plt.rc('font', family='serif')
plt.rc('xtick', labelsize='20')
plt.rc('ytick', labelsize='20')
#%%
file_dir = os.path.dirname(os.path.realpath(__file__))
idx = 2
heig = str(np.arange(5, 105, 10)[idx]).zfill(3)
files_BW = sorted(glob.glob(file_dir + "/HS/BW_flux_y%s_*.npy"%heig))
files_UW = sorted(glob.glob(file_dir + "/HS/UW*.npy"))

# BW_flux_y%s_%s.npy"%(heig, sqe)

#%%
flux_BW = np.load(files_BW[3])
flux_UW = np.load(files_UW[3])
#%%
fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(1,1,1)
ax.set_aspect(aspect=1)
con = plt.contourf(np.mean(flux_BW, axis = 0), cmap = 'gnuplot2')
# con = plt.contourf(flux_BW[1,:], levels = np.arange(-5,0.5,0.1), cmap = 'gnuplot2')
ax.set_title('', fontsize=18)
cbar = fig.colorbar(con)
cbar.ax.set_ylabel('flux',fontsize = 18)
plt.tight_layout()
# %%
fig1 = plt.figure(figsize=(8, 6))
ax = fig1.add_subplot(1,1,1)
ax.plot(flux_BW[1,:, 1024], 'k--', label='N-off2',linewidth = 1)
# labels = np.arange(0,9.6,1)
# plt.yticks(labels, labels)
ax.set_yscale('symlog')
legend = plt.legend(loc='best', frameon=False,fontsize = 18)
plt.grid(linestyle = ':')
plt.tight_layout()
# plt.savefig('output_plots/...')

# %%
