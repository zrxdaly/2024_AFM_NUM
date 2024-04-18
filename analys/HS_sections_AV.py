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

files_B = sorted(glob.glob(file_dir + "/HS/PAV_B015_*.npy"))
files_U = sorted(glob.glob(file_dir + "/HS/PAV_M015_*.npy"))
files_TV = sorted(glob.glob(file_dir + "/HS/PAV_TV015_*.npy"))
# %%
B_HS_OP = np.mean(np.load(files_B[0]), axis = 0)*273.15/9.81
B_HS_RE = np.mean(np.load(files_B[1]), axis = 0)*273.15/9.81

TV_HS_OP = np.mean(np.load(files_TV[0]), axis = 0)-273.15
TV_HS_RE = np.mean(np.load(files_TV[1]), axis = 0)-273.15

U_HS_OP = np.mean(np.load(files_U[0]), axis = 0)
U_HS_RE = np.mean(np.load(files_U[1]), axis = 0)
# %%
SA = 30
b2_OP_along = B_HS_OP[:, 1024-SA:1024+SA].mean(axis = 1)
b2_RE_along = B_HS_RE[:, 1024-SA:1024+SA].mean(axis = 1)

U2_OP_along = U_HS_OP[:, 1024-SA:1024+SA].mean(axis = 1)
U2_RE_along = U_HS_RE[:, 1024-SA:1024+SA].mean(axis = 1)

b2_OP_across = B_HS_OP[500-SA:500+SA, :].mean(axis = 0)
b2_RE_across = B_HS_RE[500-SA:500+SA, :].mean(axis = 0)

U2_OP_across = U_HS_OP[500-SA:500+SA, :].mean(axis = 0)
U2_RE_across = U_HS_RE[500-SA:500+SA, :].mean(axis = 0)
#%%
TV_HS_OP[TV_HS_RE>0.99] = np.nan
TV_HS_RE[TV_HS_RE>0.99] = np.nan
# %%
TV2_OP_along = np.nanmean(TV_HS_OP[:, 1024-SA:1024+SA], axis = 1)
TV2_RE_along = np.nanmean(TV_HS_RE[:, 1024-SA:1024+SA], axis = 1)

TV2_OP_across = np.nanmean(TV_HS_OP[500-SA:500+SA, :], axis = 0)
TV2_RE_across = np.nanmean(TV_HS_RE[500-SA:500+SA, :], axis = 0)
# %%
xb_along = np.arange(2049)-500
xb_across = np.arange(2049)-1024
fig, (ax1, ax2) = plt.subplots(2, 1, figsize = (10, 6))
ax1.plot(xb_along, b2_OP_along, 'r-', label=r'$\mathrm{T_a}$ on',linewidth = 2)
ax1.plot(xb_along, b2_RE_along, 'b-', label=r'$\mathrm{T_a}$ off',linewidth = 2)
# ax1.plot(xb_along, U2_OP_along*5, 'r+', label=r'$\mathrm{T_a}$ on',linewidth = 2)
# ax1.plot(xb_along, U2_RE_along*5, 'b+', label=r'$\mathrm{T_a}$ off',linewidth = 2)
ax1.plot(xb_along, TV2_OP_along, 'r--', label=r'$\mathrm{T_p}$ on',linewidth = 2)
ax1.plot(xb_along, TV2_RE_along, 'b--', label=r'$\mathrm{T_p}$ off',linewidth = 2)
ylabel_ax1 = np.arange(1,8,2)
xlabels_ax1 = np.arange(-500,2000,500)
plt.setp(ax1, xticks=xlabels_ax1, yticks=ylabel_ax1, xlim = (-500,1500), ylim = (0, 7))

ax1.legend(loc='best', frameon=False,fontsize = 20, ncols = 2)
ax1.vlines(x = 0, ymin = -1, ymax = 8, color = "k", linewidth = 1)
ax1.set_xlabel("x streamwise distance to the WM [m]",fontsize=20)
ax1.arrow(110-500, 6, 200, 0, width = 0.2, head_width = 0.5, head_length = 30, linestyle = "--", linewidth = 1.5, edgecolor = "k", facecolor ="None")
ax1.text(x = 140-500, y = 6.3, s = "wind", fontsize = 20, c = "k")
ax1.text(x = 45-500, y = 6, s = r"\textbf{a)}", fontsize = 20, c = "k")

ax2.plot(xb_across, b2_OP_across, 'r-', label=r'$\mathrm{T_a}$ on',linewidth = 2)
ax2.plot(xb_across, b2_RE_across, 'b-', label=r'$\mathrm{T_a}$ off',linewidth = 2)
ax2.plot(xb_across, TV2_OP_across, 'r--', label=r'$\mathrm{T_p}$ on',linewidth = 2)
ax2.plot(xb_across, TV2_RE_across, 'b--', label=r'$\mathrm{T_p}$ off',linewidth = 2)
ylabel_ax1 = np.arange(1,8,2)
xlabels_ax2 = np.arange(-1000,2000,500)
plt.setp(ax2, xticks=xlabels_ax2, yticks=ylabel_ax1, xlim = (-1000,1000), ylim = (0, 7))

ax2.vlines(x = 0, ymin = -1, ymax = 8, color = "k", linewidth = 1)
ax2.set_xlabel("y cross-stream distance to the WM [m]",fontsize=20)
ax2.plot(140-1000, 6.3, fillstyle="none", markersize=15, color='k')
ax2.plot(140-1000, 6.3, marker="x", markersize=11, color='k')
ax2.text(x = 170-1000, y = 6.1, s = "wind", fontsize = 20, c = "k")
ax2.text(x = 45-1000, y = 6, s = r"\textbf{b)}", fontsize = 20, c = "k")

fig.supylabel('temperature [$^\circ$C]',fontsize=20)
plt.tight_layout()
plt.savefig('Sections_AV.pdf')
#%%
boltz =  5.67E-8
Gconst =  9.81
T_ref =  273.15
Kd = 0.024 
VF_s =  0.1
VF_g =  VF_s
VF_l =  (1. - VF_s) 
eps_s =  0.8
eps_g =  0.98
eps_l =   0.96 
T_s = 268.15
T_g  = 272.15
def LWnet(TV):
    TV = TV + 273.15
    Lwin = 0.5 * VF_s * eps_s * boltz * T_s ** 4 + 0.5 * VF_g * eps_g * boltz * T_g ** 4 + 1. * VF_l * eps_l * boltz * TV**4
    Lwout = eps_l * boltz * TV**4
    Lwnet = Lwin - Lwout
    return Lwnet

rho_a = 1.27 
Cp_a = 1005.
L_l = 8E-2
vis = 1.718E-5 / rho_a
def Hheat(TV, TA, U):
    TV = TV + 273.15
    T_a = TA + T_ref
    gstar = Gconst * (TV - T_a) / T_a
    M = np.sqrt(U**2 + np.fabs(2 * L_l * gstar))
    Re = M * L_l / vis
    Nu = np.where(Re > 2E4, 0.032 * pow(Re, 0.8), 0.6 * pow(Re, 0.5))
    rH = L_l / Nu / Kd * Cp_a * rho_a
    H = Cp_a * rho_a / rH * (TV - T_a)
    return H
#%%
Lwnet_OP = LWnet(TV_HS_OP)
Lwnet_RE = LWnet(TV_HS_RE)

H_OP = Hheat(TV_HS_OP, B_HS_OP, U_HS_OP)
H_RE = Hheat(TV_HS_RE, B_HS_RE, U_HS_RE)

# %%
Lwnet_OP_along = np.nanmean(Lwnet_OP[:, 1024-SA:1024+SA], axis = 1)
Lwnet_RE_along = np.nanmean(Lwnet_RE[:, 1024-SA:1024+SA], axis = 1)

H_OP_along = np.nanmean(H_OP[:, 1024-SA:1024+SA], axis = 1)
H_RE_along = np.nanmean(H_RE[:, 1024-SA:1024+SA], axis = 1)

Lwnet_OP_across = np.nanmean(Lwnet_OP[500-SA:500+SA, :], axis = 0)
Lwnet_RE_across = np.nanmean(Lwnet_RE[500-SA:500+SA, :], axis = 0)

H_OP_across = np.nanmean(H_OP[500-SA:500+SA, :], axis = 0)
H_RE_across = np.nanmean(H_RE[500-SA:500+SA, :], axis = 0)

# %%
xb_along = np.arange(2049)-500
xb_across = np.arange(2049)-1024
fig, (ax1, ax2) = plt.subplots(2, 1, figsize = (11, 6))
ax1.plot(xb_along, H_OP_along, 'r-', label='$H$ on',linewidth = 2)
ax1.plot(xb_along, H_RE_along, 'b-', label='$H$ off',linewidth = 2)
ax1.plot(xb_along, Lwnet_OP_along, 'r--', label='$R_{n e t}$ on',linewidth = 2)
ax1.plot(xb_along, Lwnet_RE_along, 'b--', label='$R_{n e t}$ off',linewidth = 2)
ylabel_ax1 = np.arange(-8,-2,2)
xlabels_ax1 = np.arange(-500,2000,500)
plt.setp(ax1, xticks=xlabels_ax1, yticks=ylabel_ax1, xlim = (-500,1500), ylim = (-8,-2))

ax1.legend(loc='best', frameon=False,fontsize = 20, ncols = 2)
ax1.vlines(x = 0, ymin = -10, ymax = 0, color = "k", linewidth = 1)
ax1.set_xlabel("x streamwise distance to the WM [m]",fontsize=20)
ax1.arrow(110-500, -1.2-6, 200, 0, width = 0.2, head_width = 0.5, head_length = 30, linestyle = "--", linewidth = 1.5, edgecolor = "k", facecolor ="None")
ax1.text(x = 140-500, y = -0.9-6, s = "wind", fontsize = 20, c = "k")
ax1.text(x = 45-500, y = -1.2-6, s = r"\textbf{c)}", fontsize = 20, c = "k")

# axins1 = ax1.inset_axes([0.60, 0.60, 0.40, 0.40])
# axins1.plot(xb_along, H_OP_along, 'r-', label='$H$ on',linewidth = 2)
# axins1.plot(xb_along, H_RE_along, 'b-', label='$H$ off',linewidth = 2)
# axins1.plot(xb_along, Lwnet_OP_along, 'r--', label='Lwnet on',linewidth = 2)
# axins1.plot(xb_along, Lwnet_RE_along, 'b--', label='Lwnet off',linewidth = 2)
# axins1.tick_params(top=True, labeltop=True, bottom=False, labelbottom=False, labelsize = "small")
# # setting of a zoomed graph
# x1, x2, y1, y2 = 500, 1500, -4.0, -3.5
# axins1.set_xlim(x1, x2)
# axins1.set_ylim(y1, y2)
# ax1.indicate_inset_zoom(axins1, alpha = 1)

ax2.plot(xb_across, H_OP_across, 'r-', label='H on',linewidth = 2)
ax2.plot(xb_across, H_RE_across, 'b-', label='H off',linewidth = 2)
ax2.plot(xb_across, Lwnet_OP_across, 'r--', label='$R_{n e t}$ on',linewidth = 2)
ax2.plot(xb_across, Lwnet_RE_across, 'b--', label='$R_{n e t}$ off',linewidth = 2)
ylabel_ax1 = np.arange(-8,-2,2)
xlabels_ax2 = np.arange(-1000,2000,500)
plt.setp(ax2, xticks=xlabels_ax2, yticks=ylabel_ax1, xlim = (-1000,1000), ylim = (-8,-2))

ax2.vlines(x = 0, ymin = -10, ymax = 0, color = "k", linewidth = 1)
ax2.set_xlabel("y cross-stream distance to the WM [m]",fontsize=20)
ax2.plot(140-1000, -1.-6, fillstyle="none", markersize=15, color='k')
ax2.plot(140-1000, -1.-6, marker="x", markersize=11, color='k')
ax2.text(x = 170-1000, y = -1.2-6, s = "wind", fontsize = 20, c = "k")
ax2.text(x = 45-1000, y = -1.2-6, s = r"\textbf{d)}", fontsize = 20, c = "k")

# axins = ax2.inset_axes([0.60, 0.60, 0.40, 0.40])
# axins.plot(xb_across, H_OP_across, 'r-', label='H on',linewidth = 2)
# axins.plot(xb_across, H_RE_across, 'b-', label='H off',linewidth = 2)
# axins.plot(xb_across, Lwnet_OP_across, 'r--', label='$R_{n e t}$ on',linewidth = 2)
# axins.plot(xb_across, Lwnet_RE_across, 'b--', label='$R_{n e t}$ off',linewidth = 2)
# axins.tick_params(top=True, labeltop=True, bottom=False, labelbottom=False, labelsize = "small")
# # setting of a zoomed graph
# x1, x2, y1, y2 = 400, 1000, -4.0, -3.5
# axins.set_xlim(x1, x2)
# axins.set_ylim(y1, y2)
# ax2.indicate_inset_zoom(axins, alpha = 1)

fig.supylabel('energy budget [$W m^{-2}$]',fontsize=20)
plt.tight_layout()
plt.savefig('Sections_EB.pdf')
# %% switch the plot arrangement by plotting streamwise together and cross-stream together
xb_along = np.arange(2049)-500
xb_across = np.arange(2049)-1024
fig, (ax1, ax2) = plt.subplots(2, 1, figsize = (10, 6), sharex=True)
ax1.plot(xb_along, b2_OP_along, 'r-', label=r'$\mathrm{T_a}$ on',linewidth = 2)
ax1.plot(xb_along, b2_RE_along, 'b-', label=r'$\mathrm{T_a}$ off',linewidth = 2)
ax1.plot(xb_along, TV2_OP_along, 'r--', label=r'$\mathrm{T_p}$ on',linewidth = 2)
ax1.plot(xb_along, TV2_RE_along, 'b--', label=r'$\mathrm{T_p}$ off',linewidth = 2)
ylabel_ax1 = np.arange(1,8,2)
xlabels_ax1 = np.arange(-500,2000,500)
plt.setp(ax1, xticks=xlabels_ax1, yticks=ylabel_ax1, xlim = (-500,1500), ylim = (0, 7))

ax1.legend(loc='best', frameon=False,fontsize = 20, ncols = 2)
ax1.vlines(x = 0, ymin = -1, ymax = 8, color = "k", linewidth = 1)
# ax1.set_xlabel("x streamwise distance to the WM [m]",fontsize=20)
ax1.arrow(110-500, 6, 200, 0, width = 0.2, head_width = 0.5, head_length = 30, linestyle = "--", linewidth = 1.5, edgecolor = "k", facecolor ="None")
ax1.text(x = 140-500, y = 6.3, s = "wind", fontsize = 20, c = "k")
ax1.text(x = 45-500, y = 6, s = r"\textbf{a)}", fontsize = 20, c = "k")
ax1.set_ylabel('temperature [$^\circ$C]',fontsize=20)

ax2.plot(xb_along, H_OP_along, 'r-', label='$H$ on',linewidth = 2)
ax2.plot(xb_along, H_RE_along, 'b-', label='$H$ off',linewidth = 2)
ax2.plot(xb_along, Lwnet_OP_along, 'r--', label='$R_{n e t}$ on',linewidth = 2)
ax2.plot(xb_along, Lwnet_RE_along, 'b--', label='$R_{n e t}$ off',linewidth = 2)
ylabel_ax2 = np.arange(-8,-2,2)
xlabels_ax2 = np.arange(-500,2000,500)
plt.setp(ax2, xticks=xlabels_ax2, yticks=ylabel_ax2, xlim = (-500,1500), ylim = (-8,-2))

ax2.legend(loc='best', frameon=False,fontsize = 20, ncols = 2)
ax2.vlines(x = 0, ymin = -10, ymax = 0, color = "k", linewidth = 1)
ax2.set_xlabel("x streamwise distance to the WM [m]",fontsize=20)
ax2.arrow(110-500, -1.2-6, 200, 0, width = 0.2, head_width = 0.5, head_length = 30, linestyle = "--", linewidth = 1.5, edgecolor = "k", facecolor ="None")
ax2.text(x = 140-500, y = -0.9-6, s = "wind", fontsize = 20, c = "k")
ax2.text(x = 45-500, y = -1.2-6, s = r"\textbf{b)}", fontsize = 20, c = "k")
ax2.set_ylabel('energy budget [$W m^{-2}$]',fontsize=20)

# fig.supylabel('temperature [$^\circ$C]',fontsize=20)
plt.tight_layout()
plt.savefig('Sections_SW.pdf')
# %%
xb_along = np.arange(2049)-500
xb_across = np.arange(2049)-1024
fig, (ax1, ax2) = plt.subplots(2, 1, figsize = (10, 6), sharex=True)
ax1.plot(xb_across, b2_OP_across, 'r-', label=r'$\mathrm{T_a}$ on',linewidth = 2)
ax1.plot(xb_across, b2_RE_across, 'b-', label=r'$\mathrm{T_a}$ off',linewidth = 2)
ax1.plot(xb_across, TV2_OP_across, 'r--', label=r'$\mathrm{T_p}$ on',linewidth = 2)
ax1.plot(xb_across, TV2_RE_across, 'b--', label=r'$\mathrm{T_p}$ off',linewidth = 2)
ylabel_ax1 = np.arange(1,8,2)
xlabels_ax1 = np.arange(-1000,2000,500)
plt.setp(ax1, xticks=xlabels_ax1, yticks=ylabel_ax1, xlim = (-1000,1000), ylim = (0, 7))

ax1.legend(loc='best', frameon=False,fontsize = 20, ncols = 2)
ax1.vlines(x = 0, ymin = -1, ymax = 8, color = "k", linewidth = 1)
# ax1.set_xlabel("y cross-stream distance to the WM [m]",fontsize=20)
ax1.plot(140-1000, 6.3, fillstyle="none", markersize=15, color='k')
ax1.plot(140-1000, 6.3, marker="x", markersize=11, color='k')
ax1.text(x = 170-1000, y = 6.1, s = "wind", fontsize = 20, c = "k")
ax1.text(x = 45-1000, y = 6, s = r"\textbf{c)}", fontsize = 20, c = "k")
ax1.set_ylabel('temperature [$^\circ$C]',fontsize=20)

ax2.plot(xb_across, H_OP_across, 'r-', label='H on',linewidth = 2)
ax2.plot(xb_across, H_RE_across, 'b-', label='H off',linewidth = 2)
ax2.plot(xb_across, Lwnet_OP_across, 'r--', label='$R_{n e t}$ on',linewidth = 2)
ax2.plot(xb_across, Lwnet_RE_across, 'b--', label='$R_{n e t}$ off',linewidth = 2)
ylabel_ax1 = np.arange(-8,-2,2)
xlabels_ax2 = np.arange(-1000,2000,500)
plt.setp(ax2, xticks=xlabels_ax2, yticks=ylabel_ax1, xlim = (-1000,1000), ylim = (-8,-2))

ax2.legend(loc='best', frameon=False,fontsize = 20, ncols = 2)
ax2.vlines(x = 0, ymin = -10, ymax = 0, color = "k", linewidth = 1)
ax2.set_xlabel("y cross-stream distance to the WM [m]",fontsize=20)
ax2.plot(140-1000, -1.-6, fillstyle="none", markersize=15, color='k')
ax2.plot(140-1000, -1.-6, marker="x", markersize=11, color='k')
ax2.text(x = 170-1000, y = -1.2-6, s = "wind", fontsize = 20, c = "k")
ax2.text(x = 45-1000, y = -1.2-6, s = r"\textbf{d)}", fontsize = 20, c = "k")
ax2.set_ylabel('energy budget [$W m^{-2}$]',fontsize=20)

plt.tight_layout()
plt.savefig('Sections_CS.pdf')

# %%
