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
plt.style.use('science')
#%%
file_dir = os.path.dirname(os.path.realpath(__file__))

files_b = sorted(glob.glob(file_dir + "/TZ/TZ_b_*.npy"))
files_u = sorted(glob.glob(file_dir + "/TZ/TZ_U_*.npy"))
files_v = sorted(glob.glob(file_dir + "/TZ/TZ_V_*.npy"))
files_w = sorted(glob.glob(file_dir + "/TZ/TZ_W_*.npy"))
#%%
b_section = np.load(files_b[2])
for i in range(3, len(files_b)):
    b_section0 = np.load(files_b[i])
    b_section = np.append(b_section, b_section0, axis = 0)
    
u_section = np.load(files_u[2])
for i in range(3, len(files_u)):
    u_section0 = np.load(files_u[i])
    u_section = np.append(u_section, u_section0, axis = 0)
    
v_section = np.load(files_v[2])
for i in range(3, len(files_v)):
    v_section0 = np.load(files_v[i])
    v_section = np.append(v_section, v_section0, axis = 0)
    
w_section = np.load(files_w[2])
for i in range(3, len(files_w)):
    w_section0 = np.load(files_w[i])
    w_section = np.append(w_section, w_section0, axis = 0)
#%%
PAV_b_OP = b_section[72*0:72*(0+1), :]
PAV_u_OP = u_section[72*0:72*(0+1), :]
PAV_v_OP = v_section[72*0:72*(0+1), :]
PAV_w_OP = w_section[72*0:72*(0+1), :]
for j in range(1, 8):
    print(j)
    PAV_b_OP = PAV_b_OP + b_section[72*j:72*(j+1), :]
    PAV_u_OP = PAV_u_OP + u_section[72*j:72*(j+1), :]
    PAV_v_OP = PAV_v_OP + v_section[72*j:72*(j+1), :]
    PAV_w_OP = PAV_w_OP + w_section[72*j:72*(j+1), :]
PAV_b_OP = PAV_b_OP/8
PAV_u_OP = PAV_u_OP/8
PAV_v_OP = PAV_v_OP/8
PAV_w_OP = PAV_w_OP/8
#%%
for ind in range(72):
    with plt.ioff():
        TT = str(ind).zfill(2)
        b_05 = PAV_b_OP[ind,:]
        u_05 = PAV_u_OP[ind,:]
        v_05 = PAV_v_OP[ind,:]
        w_05 = PAV_w_OP[ind,:]
        z0 = 1024-200
        z1 = 1024+200
        H_lim = 21
        z = np.linspace(z0, z1, 401)
        y = np.linspace(0, H_lim-1, H_lim)
        Z, Y = np.meshgrid(z, y)
        fig, ax = plt.subplots(figsize=(16, 4))
        con = ax.contourf(Z, Y, b_05[:H_lim, z0:z1+1], levels = np.arange(1,11,0.1), cmap="gnuplot2")
        tpp1 = ax.contour(Z, Y, b_05[:H_lim, z0:z1+1], levels = np.arange(1,11,1), colors = 'w', alpha = 0.9, linestyles = "solid")
        ax.clabel(tpp1, fmt='%2.1f', colors='w', fontsize=12)
        for c in con.collections:
            c.set_edgecolor("face")
        # Create wind direction plot on top
        Dx_Arr = 10
        ax.quiver(Z[:,::Dx_Arr], Y[:,::Dx_Arr], u_05[:H_lim, z0:z1+1:Dx_Arr], w_05[:H_lim, z0:z1+1:Dx_Arr], angles='uv', scale = 0.3, scale_units="x", units = "x", headwidth = 2)
        # Set plot title and labels 
        # ax.set_title('Temperature Contour Plot with Wind Direction')
        # ax.set_xlabel('distance [m]')
        # ax.set_ylabel('height [m]')
        plt.tight_layout()
        plt.savefig(file_dir + "/TZ/movies_PTA/TZ_t%s.png"%TT)
        # break
        # plt.close()
#%    
# %%
