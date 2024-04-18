#%% 
import numpy as np
import glob
import os
import matplotlib.pyplot as plt
plt.rc('font', family='serif')
plt.rc('xtick', labelsize='20')
plt.rc('ytick', labelsize='20')
plt.rc('text', usetex=True)
import scienceplots
plt.style.use(['science', 'grid', 'scatter'])

from matplotlib import cm
import matplotlib.colors as colors
import matplotlib.colors as mcolors
cmap = cm.coolwarm_r
#%%
file_dir = os.path.dirname(os.path.realpath(__file__))
Rb_file = sorted(glob.glob(file_dir + "/Radial_profile/Rb_AG_OP.npy"))
RM_file = sorted(glob.glob(file_dir + "/Radial_profile/RM_AG_OP.npy"))
#%%
Rb_AG= np.load(Rb_file[0])[:,:250]
RM_AG = np.load(RM_file[0])[:,:250]
# %% plotting together
def plotting(Rb_AG1_T1, RM_AG1_T1, zone):
    Rb_AG1_T1_F = Rb_AG1_T1.flatten()
    RM_AG1_T1_F = RM_AG1_T1.flatten()
    R_l = Rb_AG1_T1.shape[1]
    R_range = np.arange(1, R_l*2, 2)
    liftH = 0.2
    Heig_N = np.arange(liftH, 21 + liftH, 1.)
    Heig_C = np.arange(0.5, 21, 1)

    dis_s_AG1= np.tile(R_range, 21)
    Heig_C =  np.arange(0.5, 21, 1)
    H_s_AG1 = np.repeat(Heig_C, R_l)
    liftH = 0.2
    d_rate = liftH/10
    for i in range(11):
        H_s_AG1[i::R_l] = H_s_AG1[i::R_l] + liftH - d_rate * i
    s_norm = mcolors.TwoSlopeNorm(vmin=0., vcenter=20., vmax=R_range[-1])
    ln = 1.5
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize = (10, 8), sharey=True)
    im = ax1.scatter(RM_AG1_T1_F, H_s_AG1, (1.2-s_norm(dis_s_AG1))*30, dis_s_AG1, norm=colors.TwoSlopeNorm(vcenter=10), cmap = cmap, marker="o")
    ax1.plot(RM_AG1_T1[:, 0],  Heig_C, 'r--',linewidth = ln, label = "1 m", alpha = 0.8)
    ax1.plot(RM_AG1_T1[:, 10],  Heig_C, 'k--',linewidth = ln, label = "21 m", alpha = 0.8)
    ax1.plot(RM_AG1_T1[:, -1],  Heig_C, 'b--',linewidth = ln, label = "%s m"%(str(R_range[-1])), alpha = 1)
    ax1.set_xlabel("$M$ [$m~s^{-1}$]",fontsize=18)
    ax1.set_ylabel('Height [m]',fontsize=18)
    ylabels = np.arange(0,20,2)
    xlabels1 = np.arange(-0,12,2)
    plt.setp(ax1, xticks=xlabels1, yticks=ylabels, xlim = (-0,12), ylim = (0, 20))

    axins = ax1.inset_axes([0.65, 0.04, 0.25, 0.35])
    # axins.plot(RM_AG1_T1[:, 0],  Heig_N, 'r--',linewidth = ln, label = "1 m", alpha = 0.8)
    axins.plot(RM_AG1_T1[:, 10],  Heig_C, 'k--',linewidth = ln, label = "21 m", alpha = 0.8)
    axins.plot(RM_AG1_T1[:, -1],  Heig_C, 'b--',linewidth = ln, label = "%s m"%(str(R_range[-1])), alpha = 1)
    # setting of a zoomed graph
    x1, x2, y1, y2 = 0, 2, 0, 20
    axins.set_xlim(x1, x2)
    axins.set_ylim(y1, y2)
    ax1.indicate_inset_zoom(axins)

    im = ax2.scatter(Rb_AG1_T1_F, H_s_AG1, (1.2-s_norm(dis_s_AG1))*30, dis_s_AG1, norm=colors.TwoSlopeNorm(vcenter=10), cmap = cmap, marker="o")
    ax2.plot(Rb_AG1_T1[:, 0],  Heig_C, 'r--',linewidth = ln, label = "1 m", alpha = 0.8)
    ax2.plot(Rb_AG1_T1[:, 10],  Heig_C, 'k--',linewidth = ln, label = "21 m", alpha = 0.8)
    ax2.plot(Rb_AG1_T1[:, -1],  Heig_C, 'b--',linewidth = ln, label = "%s m"%(str(R_range[-1])), alpha = 1)
    # ax2.set_ylabel('Height [m]',fontsize=18)
    ax2.set_xlabel("$T_a$ [$^\circ$C]",fontsize=18)
    ax2.legend(loc='best', frameon=False,fontsize = 20)
    xlabels2 = np.arange(0,10,2)
    plt.setp(ax2, xticks=xlabels2, yticks=ylabels, xlim = (0,10), ylim = (0, 20))
    fig.subplots_adjust(bottom=0.1, top=0.9, left=0.1, right=0.8,
                        wspace=0.02, hspace=0.05)
    cb_ax = fig.add_axes([0.83, 0.1, 0.02, 0.8])
    cbar = fig.colorbar(im, cax=cb_ax)
    cbar.ax.set_title('r [m]', fontsize = 20)
    cbar.set_ticks([1, 5, 20, 100, 200, 400])
    ax1.text(x = 10, y = 19, s = r"\textbf{a)}", fontsize = 20, c = "k")
    ax2.text(x = 8,  y = 19, s = r"\textbf{b)}", fontsize = 20, c = "k")
    # cb = fig1.colorbar(im, shrink=0.8,aspect=60, fraction=0.046, pad=0.04)
    # cb.set_label('Distance from WM')
    # cb.ax.set_title('r[m]')
    # cb.ax1.set_xlabel('r [m]', fontsize = 20)
    # plt.tight_layout()
    
    plt.savefig(file_dir + "/Radial_profile/%s.pdf"%zone, bbox_inches='tight')
# %%
plotting(Rb_AG, RM_AG, "R_Mb")
#%%







# %% plot the difference 
Rb_file_RE = sorted(glob.glob(file_dir + "/radial/radial_all/Rb_AG_RE.npy"))
RM_file_RE = sorted(glob.glob(file_dir + "/radial/radial_all/RM_AG_RE.npy"))
#%%
Rb_AG_RE = np.load(Rb_file_RE[0])[:,:250]
RM_AG_RE = np.load(RM_file_RE[0])[:,:250]
# %% plotting difference
Rb_DT = Rb_AG - Rb_AG_RE
RM_DM = RM_AG - RM_AG_RE
#%%
Rb_DT = Rb_AG - Rb_AG[:,-1][:, np.newaxis]
RM_DM = RM_AG - RM_AG[:,-1][:, np.newaxis]
#%%
plt.contourf(RM_DM, levels = np.arange(-0.1, 5, 0.01), cmap = "gnuplot2")
# %% plotting together
def plotting_DF(Rb_AG1_T1, RM_AG1_T1, zone):
    Rb_AG1_T1_F = Rb_AG1_T1.flatten()
    RM_AG1_T1_F = RM_AG1_T1.flatten()
    R_l = Rb_AG1_T1.shape[1]
    R_range = np.arange(1, R_l*2, 2)
    Heig_C = np.arange(0, 21, 1)

    dis_s_AG1= np.tile(R_range, 21)
    Heig_C = np.arange(0., 21., 1.)
    H_s_AG1 = np.repeat(Heig_C, R_l)
    s_norm = mcolors.TwoSlopeNorm(vmin=0., vcenter=20., vmax=R_range[-1])
    ln = 1.5
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize = (10, 8), sharey=True)
    ax1.plot(RM_AG1_T1[:, 1],  Heig_C, 'r--',linewidth = ln, label = "0 m", alpha = 0.8)
    ax1.plot(RM_AG1_T1[:, 6],  Heig_C, 'r--',linewidth = ln, label = "13 m", alpha = 0.8)
    ax1.plot(RM_AG1_T1[:, 8],  Heig_C, 'r--',linewidth = ln, label = "17 m", alpha = 0.8)
    ax1.plot(RM_AG1_T1[:, 10],  Heig_C, 'k--',linewidth = ln, label = "21 m", alpha = 0.8)
    ax1.plot(RM_AG1_T1[:, 25],  Heig_C, 'k--',linewidth = ln, label = "21 m", alpha = 0.8)
    ax1.plot(RM_AG1_T1[:, 50],  Heig_C, 'k--',linewidth = ln, label = "101 m", alpha = 0.8)
    ax1.plot(RM_AG1_T1[:, 75],  Heig_C, 'k--',linewidth = ln, label = "151 m", alpha = 0.8)
    ax1.plot(RM_AG1_T1[:, 100],  Heig_C, 'k--',linewidth = ln, label = "201 m", alpha = 0.8)
    ax1.plot(RM_AG1_T1[:, 125],  Heig_C, 'k--',linewidth = ln, label = "251 m", alpha = 0.8)
    ax1.plot(RM_AG1_T1[:, 150],  Heig_C, 'k--',linewidth = ln, label = "301 m", alpha = 0.8)
    ax1.plot(RM_AG1_T1[:, -1],  Heig_C, 'b--',linewidth = ln, label = "%s m"%(str(R_range[-1])), alpha = 1)
    ax1.set_xlabel("M [m s$^{-1}$]",fontsize=18)
    ax1.set_ylabel('Height [m]',fontsize=18)
    ylabels = np.arange(0,21,2)
    # xlabels1 = np.arange(-2,10,2)
    # plt.setp(ax1, xticks=xlabels1, yticks=ylabels, xlim = (-2,10), ylim = (0, 20))
    plt.setp(ax1, yticks=ylabels, ylim = (0, 20))

    # im = ax2.scatter(Rb_AG1_T1_F, H_s_AG1, (1.2-s_norm(dis_s_AG1))*30, dis_s_AG1, norm=colors.TwoSlopeNorm(vcenter=10), cmap = cmap, marker="o")
    ax2.plot(Rb_AG1_T1[:, 3],  Heig_C, 'r--',linewidth = ln, label = "1 m", alpha = 0.8)
    ax2.plot(Rb_AG1_T1[:, 6],  Heig_C, 'r--',linewidth = ln, label = "13 m", alpha = 0.8)
    ax2.plot(Rb_AG1_T1[:, 8],  Heig_C, 'r--',linewidth = ln, label = "17 m", alpha = 0.8)
    ax2.plot(Rb_AG1_T1[:, 10],  Heig_C, 'k--',linewidth = ln, label = "21 m", alpha = 0.8)
    ax2.plot(Rb_AG1_T1[:, 25],  Heig_C, 'k--',linewidth = ln, label = "21 m", alpha = 0.8)
    ax2.plot(Rb_AG1_T1[:, 50],  Heig_C, 'k--',linewidth = ln, label = "101 m", alpha = 0.8)
    ax2.plot(Rb_AG1_T1[:, 75],  Heig_C, 'k--',linewidth = ln, label = "151 m", alpha = 0.8)
    ax2.plot(Rb_AG1_T1[:, 100],  Heig_C, 'k--',linewidth = ln, label = "201 m", alpha = 0.8)
    ax2.plot(Rb_AG1_T1[:, 125],  Heig_C, 'k--',linewidth = ln, label = "251 m", alpha = 0.8)
    ax2.plot(Rb_AG1_T1[:, 150],  Heig_C, 'k--',linewidth = ln, label = "301 m", alpha = 0.8)
    ax2.plot(Rb_AG1_T1[:, -1],  Heig_C, 'b--',linewidth = ln, label = "%s m"%(str(R_range[-1])), alpha = 1)
    # ax2.set_ylabel('Height [m]',fontsize=18)
    ax2.set_xlabel("Ta [$^\circ$C]",fontsize=18)
    ax2.legend(loc='best', frameon=False,fontsize = 20)
    xlabels2 = np.arange(-4,6,2)
    plt.setp(ax2, xticks=xlabels2, yticks=ylabels, xlim = (-4,6), ylim = (0, 20))

    # cb = fig1.colorbar(im, shrink=0.8,aspect=60, fraction=0.046, pad=0.04)
    # cb.set_label('Distance from WM')
    # cb.ax.set_title('r[m]')
    # cb.ax1.set_xlabel('r [m]', fontsize = 20)
    # plt.tight_layout()
    # plt.savefig("radial/radial_all/%s.png"%zone, bbox_inches='tight')

#%%
plotting_DF(Rb_DT, RM_DM, "AG_DF")

# %%
