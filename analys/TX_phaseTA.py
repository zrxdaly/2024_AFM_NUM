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

# files_b = sorted(glob.glob(file_dir + "/TX/TX_b_*.npy"))
# files_u = sorted(glob.glob(file_dir + "/TX/TX_U_*.npy"))
# files_v = sorted(glob.glob(file_dir + "/TX/TX_V_*.npy"))
# files_w = sorted(glob.glob(file_dir + "/TX/TX_W_*.npy"))
# #%%
# b_section = np.load(files_b[2])
# for i in range(3, len(files_b)):
#     b_section0 = np.load(files_b[i])
#     b_section = np.append(b_section, b_section0, axis = 0)
    
# u_section = np.load(files_u[2])
# for i in range(3, len(files_u)):
#     u_section0 = np.load(files_u[i])
#     u_section = np.append(u_section, u_section0, axis = 0)
    
# w_section = np.load(files_w[2])
# for i in range(3, len(files_w)):
#     w_section0 = np.load(files_w[i])
#     w_section = np.append(w_section, w_section0, axis = 0)
# #%%
# PAV_b_OP = b_section[72*0:72*(0+1), :]
# PAV_u_OP = u_section[72*0:72*(0+1), :]
# PAV_w_OP = w_section[72*0:72*(0+1), :]
# for j in range(1, 8):
#     print(j)
#     PAV_b_OP = PAV_b_OP + b_section[72*j:72*(j+1), :]
#     PAV_u_OP = PAV_u_OP + u_section[72*j:72*(j+1), :]
#     PAV_w_OP = PAV_w_OP + w_section[72*j:72*(j+1), :]
# PAV_b_OP = PAV_b_OP/8
# PAV_u_OP = PAV_u_OP/8
# PAV_w_OP = PAV_w_OP/8
#%%
files_b = sorted(glob.glob(file_dir + "/TX/PAV_B_OP.npy"))
files_u = sorted(glob.glob(file_dir + "/TX/PAV_U_OP.npy"))
files_v = sorted(glob.glob(file_dir + "/TX/PAV_V_OP.npy"))
files_w = sorted(glob.glob(file_dir + "/TX/PAV_W_OP.npy"))

PAV_b_OP = np.load(files_b[0])*273.15/9.81
PAV_u_OP = np.load(files_u[0])
PAV_v_OP = np.load(files_v[0])
PAV_w_OP = np.load(files_w[0])

#%%
for ind in range(72):
    with plt.ioff():
        TT = str(ind).zfill(2)
        b_05 = PAV_b_OP[ind,:]
        u_05 = PAV_u_OP[ind,:]
        w_05 = PAV_w_OP[ind,:]
        H_lim = 25
        x = np.linspace(350, 750, 401)
        y = np.linspace(0, H_lim-1, H_lim)
        X, Y = np.meshgrid(x, y)
        fig, ax = plt.subplots(figsize=(16, 4))
        con = ax.contourf(X, Y, b_05.T[:H_lim, 350:751], levels = np.arange(1,11,0.1), cmap="gnuplot2")
        tpp1 = ax.contour(X, Y, b_05.T[:H_lim, 350:751], levels = np.arange(1,11,1), colors = 'w', alpha = 0.9, linestyles = "solid")
        ax.clabel(tpp1, fmt='%2.1f', colors='w', fontsize=12)
        for c in con.collections:
            c.set_edgecolor("face")
        # Create wind direction plot on top
        Dx_Arr = 15
        ax.quiver(X[:,::Dx_Arr], Y[:,::Dx_Arr], u_05.T[:H_lim, 350:751:Dx_Arr], w_05.T[:H_lim, 350:751:Dx_Arr], angles='uv', scale = 0.3, scale_units="x", units = "x", headwidth = 2)
        # Set plot title and labels 
        # ax.set_title('Temperature Contour Plot with Wind Direction')
        # ax.set_xlabel('distance [m]')
        # ax.set_ylabel('height [m]')
        plt.tight_layout()
        plt.savefig(file_dir + "/TX/movies_PAV/TX_t%s.png"%TT)
        # break
        # plt.close()
#%    
# %% the right format
H_lim = 25
x = np.linspace(350, 750, 401) -500
y = np.linspace(0, H_lim-1, H_lim)
X, Y = np.meshgrid(x, y)
fig, (ax1, ax2) = plt.subplots(2, 1, figsize = (10, 6), sharex=True)

con1 = ax1.contourf(X, Y, PAV_b_OP[40,:].T[:H_lim, 350:751], levels = np.arange(1,11,0.1), cmap="plasma")
con2 = ax2.contourf(X, Y,  PAV_b_OP[6,:].T[:H_lim, 350:751], levels = np.arange(1,11,0.1), cmap="plasma")

tpp1 = ax1.contour(X, Y, PAV_b_OP[40,:].T[:H_lim, 350:751], levels = np.arange(1,11,1), colors = 'w', alpha = 0.9, linestyles = "solid")
tpp2 = ax2.contour(X, Y,  PAV_b_OP[6,:].T[:H_lim, 350:751], levels = np.arange(1,11,1), colors = 'w', alpha = 0.9, linestyles = "solid")

ax1.clabel(tpp1, fmt='%2.1f', colors='w', fontsize=12)
ax2.clabel(tpp2, fmt='%2.1f', colors='w', fontsize=12)

for c in con1.collections:
    c.set_edgecolor("face")
for c in con2.collections:
    c.set_edgecolor("face")
    
# Create wind direction plot on top
Dx_Arr = 15
ax1.quiver(X[:,::Dx_Arr], Y[:,::Dx_Arr], PAV_u_OP[40,:].T[:H_lim, 350:751:Dx_Arr], PAV_w_OP[40,:].T[:H_lim, 350:751:Dx_Arr], angles='uv', scale = 0.3, scale_units="x", units = "x", headwidth = 2)
ax2.quiver(X[:,::Dx_Arr], Y[:,::Dx_Arr],  PAV_u_OP[6,:].T[:H_lim, 350:751:Dx_Arr],  PAV_w_OP[6,:].T[:H_lim, 350:751:Dx_Arr], angles='uv', scale = 0.3, scale_units="x", units = "x", headwidth = 2)

ax1.set_title(r"background wind $\Longrightarrow$",fontsize=20, loc= "left")
ax1.set_ylabel("z [m]",fontsize=20)
ax2.set_ylabel("z [m]",fontsize=20)
ax2.set_xlabel("x streamwise distance to the WM [m]",fontsize=20)
ax1.text(x = 731-500, y = 17.5+4, s = r"\textbf{a)}", fontsize = 20, c = "k")
ax2.text(x = 731-500, y = 17.5+4, s = r"\textbf{b)}", fontsize = 20, c = "k")

ax1.plot(500-500,11-0.5, "k*", markersize=12)
ax2.plot(500-500,11-0.5, "k*", markersize=12)

plt.tight_layout()
plt.savefig(file_dir + "/TX_sections.pdf")

# %%
