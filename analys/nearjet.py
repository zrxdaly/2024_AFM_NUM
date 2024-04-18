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
# %%
file_dir = os.path.dirname(os.path.realpath(__file__))
files_b_OP = sorted(glob.glob(file_dir + "/TX/TX_wavecheckh3.npy"))
files_U_OP = sorted(glob.glob(file_dir + "/TX/TX_U_checkh3.npy"))
files_V_OP = sorted(glob.glob(file_dir + "/TX/TX_V_checkh3.npy"))
files_W_OP = sorted(glob.glob(file_dir + "/TX/TX_W_checkh3.npy"))
# %%
b03 = np.load(files_b_OP[0])*273.15/9.81
U03 = np.load(files_U_OP[0])
V03 = np.load(files_V_OP[0])
W03 = np.load(files_W_OP[0])
M03 = np.sqrt(U03**2 + V03**2 + W03**2)

# %%
Time = np.arange(0, 3000, 4)/60
fig, (ax1, ax2) = plt.subplots(2, 1, figsize = (12, 8), sharex=True)
ax1.plot(Time, b03[:,250]-3, '-', linewidth = 2, label = r"$\overline{T_{on}}$")
ax1.plot(Time, b03[:,480]-2, '-', linewidth = 2, label = r"$\overline{T_{off}}$")
ax1.plot(Time, b03[:,520]+0, '-', linewidth = 2, label = r"$\overline{T_{off}}$")
ax1.plot(Time, b03[:,750]+5, '-', linewidth = 2, label = r"$\overline{T_{off}}$")

ax2.plot(Time, M03[:,250]-3, '-', linewidth = 2, label = r"$\overline{T_{on}}$")
ax2.plot(Time, M03[:,480]-2, '-', linewidth = 2, label = r"$\overline{T_{off}}$")
ax2.plot(Time, M03[:,520]+0, '-', linewidth = 2, label = r"$\overline{T_{off}}$")
ax2.plot(Time, M03[:,750]+3, '-', linewidth = 2, label = r"$\overline{T_{off}}$")
# labels = np.arange(0,9.6,1)
# plt.yticks(labels, labels)
# legend = plt.legend(loc='best', frameon=False,fontsize = 18)
plt.grid(linestyle = ':')
plt.tight_layout()
# plt.savefig('output_plots/...')
# %%
# %%
# T_TX_OP = np.mean(b02[150:,:], axis = 0)
# T_TX_RE = np.mean(b02[:150,:], axis = 0)
# y = np.arange(0, 2049, 1)
# x = np.arange(0, 3000, 4)/60
# X, Y = np.meshgrid(x, y)
# fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 6), sharey=True, gridspec_kw={'width_ratios': [3,1]})
# con = ax1.contourf(X, Y, b02.T, levels = np.arange(1,8.5,0.1), cmap = 'plasma')
# # conf = ax1.contour(X, Y, T_TX.T, levels = np.arange(2, 8.5, 4), colors = 'w')
# # ax1.clabel(conf, fmt='%2.1f', colors='w', fontsize = 12)
# for c in con.collections:
#     c.set_edgecolor("face")
# ax1.axvline(x = 10, c = 'w', ls = '--', lw = 1)
# ax1.set_xlabel("Time [min]", fontsize = 20)
# ax1.set_ylabel("x [m]", fontsize = 20)

# ax1.arrow(2, 1400, 0, 200, width = 0.5, head_width = 2, head_length = 50, linestyle = "--", linewidth = 1.5, edgecolor = "w", facecolor ="None")
# ax1.text(x = 2.5, y = 1500, s = "wind", fontsize = 20, c = "w")

# ax1.plot(1,500, "w*", markersize=14)

# liner = ax2.plot(T_TX_OP, y, 'r-', label = r"$\overline{T_{on}}$")
# lineb = ax2.plot(T_TX_RE, y, 'b-', label = r"$\overline{T_{off}}$")
# lineb = ax2.plot(b02[400,:], y, 'r--', label = r"$T_{ins}$")
# ax2.set_xlabel(r"$\mathrm{T_{a}}$ [$\mathrm{^\circ C}$]", fontsize = 20)
# # ax2.grid(linestyle = ":")
# ax1.text(x = 3.5, y = 1880, s = r"\textbf{off}",fontweight="bold", fontsize = 20, c = "w")
# ax1.text(x = 28,y = 1880, s = r"\textbf{on}",fontweight="bold", fontsize = 20, c = "w")
# ax1.text(x = 1, y = 50, s = r"\textbf{a)}",fontweight="bold", fontsize = 20, c = "w")
# ax2.text(x = 7, y = 50, s = r"\textbf{b)}",fontweight="bold", fontsize = 20, c = "k")
# plt.setp(ax2, xticks = np.arange(2,10,2), yticks = np.arange(0,2048,500), ylim = (0,1201))
# plt.setp(ax1, yticks = np.arange(0,2048,500))
# legend = ax2.legend(loc='upper right', frameon=False, fontsize = 16)
# cbar = fig.colorbar(con, ax=ax2)
# cbar.ax.set_ylabel(r'$\mathrm{T_{a}}$ [$\mathrm{^\circ C}$]',fontsize = 20)
# plt.tight_layout()
# plt.savefig("wavecheck.png", dpi = 300)