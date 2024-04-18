#%% this script we compare the b2 with sim and exp
import numpy as np
from glob import glob
import os
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pandas as pd
import xarray as xr

# plt.rc('font', family='serif')
plt.rc('xtick', labelsize='20')
plt.rc('ytick', labelsize='20')
plt.rc('text', usetex=True)
import scienceplots
plt.style.use('science')

# %% dir for different cases
file_dir = os.path.dirname(os.path.realpath(__file__))
dirs_OP = sorted(glob(file_dir + "/W12/*2.npy"))

# %% get b2 in 3D file
b2_OP = np.load(dirs_OP[0])

#%%
T2_OP = b2_OP[150:,:]
#%%
index_PAD = 12
Rp = 75
PAV_b_OP = T2_OP[Rp*0:Rp*(0+1), :]
for j in range(1, 8):
    print(j)
    PAV_b_OP = PAV_b_OP + T2_OP[Rp*j:Rp*(j+1), :]
PAV_b_OP = PAV_b_OP/8
#%%
b_t06 = PAV_b_OP[index_PAD,5:276,5:256]
b_t24 = PAV_b_OP[index_PAD+18*1,5:276,5:256]
b_t42 = PAV_b_OP[index_PAD+18*2,5:276,5:256]
b_t60 = PAV_b_OP[index_PAD+18*3,5:276,5:256]
#%% here we import experimental data
in_folder = '/home/dai/Documents/Reserach/Fruit/Krabbendijke_experiment_2021_0507/DTS/calibration/'
b2_T_roll = xr.open_dataset(in_folder + "/2_contour_new/T_EF_save/" + "b2_T_roll.nc", engine = "netcdf4").rename({'__xarray_dataarray_variable__': 'tmpw'})

#%% interpolatoin
OP_T_section = slice("2021-05-07T23:40:00", "2021-05-08T00:20:00")
RE_T_section = slice("2021-05-07T22:20:00", "2021-05-07T23:00:00")
b2_T_roll_OP  = b2_T_roll.sel(time =  OP_T_section)
b2_T_roll_RE  = b2_T_roll.sel(time =  RE_T_section)

# DB2 = b2_T_roll_OP.tmpw.values - b2_T_roll_RE.tmpw.values
DB2 = b2_T_roll_OP.tmpw.values
time = b2_T_roll_OP.time.values
x = b2_T_roll_OP.x.values
y = b2_T_roll_OP.y.values
DB2_inp = xr.DataArray(DB2, coords=[time, x, y], dims=["time", "x", "y"])
DB2_inp = DB2_inp.resample(time="1S").interpolate("linear").dropna(dim='time')
#%%
PAV_DB2 = DB2_inp.values[(280 + Rp * 4 * 0):(280 + Rp * 4 * 1), :]
for i in range(1, 7):
    PAV_DB2 = PAV_DB2 + DB2_inp.values[(280 + Rp * 4 * i):(280 + Rp * 4 * (i+1)), :]
PAV_DB2 = PAV_DB2/7
#%%
Dx_Arr = 10
z0 = np.linspace(0, 270, 271)
x0 = np.linspace(0, 250, 251)
X0, Z0 = np.meshgrid(x0, z0)
# , layout="constrained"
fig, ((ax1, ax2, ax3, ax4), (ax5, ax6, ax7, ax8))= plt.subplots(2, 4, figsize = (12.5, 6), sharex=True, sharey=True)
y1 = np.linspace(0,270,522)
x1 = np.linspace(0, 250, 30)
X1, Y1 = np.meshgrid(x1, y1) 

con1 = ax1.contourf(X1, Y1, PAV_DB2[72 * 0, :], levels = np.arange(1, 9, 0.1), extend = "max", cmap = "plasma")
con2 = ax2.contourf(X1, Y1, PAV_DB2[72 * 1, :], levels = np.arange(1, 9, 0.1), extend = "max", cmap = "plasma")
con3 = ax3.contourf(X1, Y1, PAV_DB2[72 * 2, :], levels = np.arange(1, 9, 0.1), extend = "max", cmap = "plasma")
con4 = ax4.contourf(X1, Y1, PAV_DB2[72 * 3, :], levels = np.arange(1, 9, 0.1), extend = "max", cmap = "plasma")

tpp1 = ax1.contour(X1, Y1, PAV_DB2[72 * 0, :], levels = np.arange(3.5, 9, 2), colors = 'k', alpha = 0.95, linestyles = "solid")
tpp2 = ax2.contour(X1, Y1, PAV_DB2[72 * 1, :], levels = np.arange(3.5, 9, 2), colors = 'k', alpha = 0.95, linestyles = "solid")
tpp3 = ax3.contour(X1, Y1, PAV_DB2[72 * 2, :], levels = np.arange(3.5, 9, 2), colors = 'k', alpha = 0.95, linestyles = "solid")
tpp4 = ax4.contour(X1, Y1, PAV_DB2[72 * 3, :], levels = np.arange(3.5, 9, 2), colors = 'k', alpha = 0.95, linestyles = "solid")

ax1.clabel(tpp1, fmt='%2.1f', colors='w', fontsize=14)
ax2.clabel(tpp2, fmt='%2.1f', colors='w', fontsize=14)
ax3.clabel(tpp3, fmt='%2.1f', colors='w', fontsize=14)
ax4.clabel(tpp4, fmt='%2.1f', colors='w', fontsize=14)

for c in con1.collections:
    c.set_edgecolor("face")
for c in con2.collections:
    c.set_edgecolor("face")
for c in con3.collections:
    c.set_edgecolor("face")
for c in con4.collections:
    c.set_edgecolor("face")
    
ax1.arrow(200, 60, -20, 33.33, width = 2, head_width = 8, linestyle = "--", edgecolor = "w", facecolor ="None")
ax2.arrow(200, 60, -20, 33.33, width = 2, head_width = 8, linestyle = "--", edgecolor = "w", facecolor ="None")
ax3.arrow(200, 60, -20, 33.33, width = 2, head_width = 8, linestyle = "--", edgecolor = "w", facecolor ="None")
ax4.arrow(200, 60, -20, 33.33, width = 2, head_width = 8, linestyle = "--", edgecolor = "w", facecolor ="None")

ax1.text(x = 170, y = 35, s = "wind", fontsize = 20, fontweight = "500", c = "w")
ax1.text(x = 150, y = 185, s = "jet", fontsize = 20, fontweight = "500", c = "w")

ax1.arrow(150, 160, -30, 50, width = 4, head_width = 10, linestyle = "-", edgecolor = "w", facecolor ="None")
ax2.arrow(150, 160,  50, 30, width = 4, head_width = 10, linestyle = "-", edgecolor = "w", facecolor ="None")
ax3.arrow(150, 160,  30,-50, width = 4, head_width = 10, linestyle = "-", edgecolor = "w", facecolor ="None")
ax4.arrow(150, 160, -50,-30, width = 4, head_width = 10, linestyle = "-", edgecolor = "w", facecolor ="None")

con5 = ax5.contourf(X0, Z0, b_t06, levels = np.arange(1,9,0.1), cmap="plasma")
con6 = ax6.contourf(X0, Z0, b_t24, levels = np.arange(1,9,0.1), cmap="plasma")
con7 = ax7.contourf(X0, Z0, b_t42, levels = np.arange(1,9,0.1), cmap="plasma")
con8 = ax8.contourf(X0, Z0, b_t60, levels = np.arange(1,9,0.1), cmap="plasma")

tpp5 = ax5.contour(X0, Z0, b_t06, levels = np.arange(3.5, 9, 2), colors = 'k', alpha = 0.95, linestyles = "solid")
tpp6 = ax6.contour(X0, Z0, b_t24, levels = np.arange(3.5, 9, 2), colors = 'k', alpha = 0.95, linestyles = "solid")
tpp7 = ax7.contour(X0, Z0, b_t42, levels = np.arange(3.5, 9, 2), colors = 'k', alpha = 0.95, linestyles = "solid")
tpp8 = ax8.contour(X0, Z0, b_t60, levels = np.arange(3.5, 9, 2), colors = 'k', alpha = 0.95, linestyles = "solid")

ax5.clabel(tpp5, fmt='%2.1f', colors='w', fontsize=12)
ax6.clabel(tpp6, fmt='%2.1f', colors='w', fontsize=12)
ax7.clabel(tpp7, fmt='%2.1f', colors='w', fontsize=12)
ax8.clabel(tpp8, fmt='%2.1f', colors='w', fontsize=12)

for c in con5.collections:
    c.set_edgecolor("face")
for c in con6.collections:
    c.set_edgecolor("face")
for c in con7.collections:
    c.set_edgecolor("face")
for c in con8.collections:
    c.set_edgecolor("face")
    
# ax5.quiver(X0[::Dx_Arr,::Dx_Arr], Z0[::Dx_Arr,::Dx_Arr], v_t06[::Dx_Arr, ::Dx_Arr], u_t06[::Dx_Arr, ::Dx_Arr], angles='uv', scale = 0.15, scale_units="x", units = "x", headwidth = 3, color = "w")
# ax6.quiver(X0[::Dx_Arr,::Dx_Arr], Z0[::Dx_Arr,::Dx_Arr], v_t24[::Dx_Arr, ::Dx_Arr], u_t24[::Dx_Arr, ::Dx_Arr], angles='uv', scale = 0.15, scale_units="x", units = "x", headwidth = 3, color = "w")
# ax7.quiver(X0[::Dx_Arr,::Dx_Arr], Z0[::Dx_Arr,::Dx_Arr], v_t42[::Dx_Arr, ::Dx_Arr], u_t42[::Dx_Arr, ::Dx_Arr], angles='uv', scale = 0.15, scale_units="x", units = "x", headwidth = 3, color = "w")
# ax8.quiver(X0[::Dx_Arr,::Dx_Arr], Z0[::Dx_Arr,::Dx_Arr], v_t60[::Dx_Arr, ::Dx_Arr], u_t60[::Dx_Arr, ::Dx_Arr], angles='uv', scale = 0.15, scale_units="x", units = "x", headwidth = 3, color = "w")

fig.text(0.5, 0.03, 'distance [m]', ha='center', fontsize=18)
fig.text(0.03, 0.5, 'distance [m]', va='center', rotation='vertical', fontsize=18)

# divider = make_axes_locatable((ax4, ax8))
# cax = divider.append_axes("right", size="5%", pad=0.05)
# hh = plt.colorbar(con1, cax=cax)
# hh.ax.set_ylabel(r'$\Delta T$ [°C]',fontsize = 18)

fig.subplots_adjust(bottom=0.1, top=0.9, left=0.1, right=0.8,
                    wspace=0.02, hspace=0.02)
cb_ax = fig.add_axes([0.83, 0.1, 0.02, 0.8])
cbar = plt.colorbar(con1, cax=cb_ax)
cbar.ax.set_ylabel(r'$\Delta T$ [°C]',fontsize = 20)

ax5.plot(150,150, "k*", markersize=8)
ax6.plot(150,150, "k*", markersize=8)
ax7.plot(150,150, "k*", markersize=8)
ax8.plot(150,150, "k*", markersize=8)

# ax.set_xlabel('distance [m]')
# ax.set_ylabel('height [m]')
plt.savefig("W12/PAV_DB2.pdf", bbox_inches='tight')

# %%
