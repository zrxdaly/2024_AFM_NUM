# %%
import numpy as np
import glob
import os 

#%%
file_dir = os.path.dirname(os.path.realpath(__file__))

files_b = sorted(glob.glob(file_dir + "/TX/TX_B_*.npy"))
files_u = sorted(glob.glob(file_dir + "/TX/TX_U_*.npy"))
files_v = sorted(glob.glob(file_dir + "/TX/TX_V_*.npy"))
files_w = sorted(glob.glob(file_dir + "/TX/TX_W_*.npy"))
# %%
b_section = np.load(files_b[6])
for i in range(7, len(files_b)):
    b_section0 = np.load(files_b[i])
    b_section = np.append(b_section, b_section0, axis = 0)
    
u_section = np.load(files_u[6])
for i in range(7, len(files_u)):
    u_section0 = np.load(files_u[i])
    u_section = np.append(u_section, u_section0, axis = 0)
    
w_section = np.load(files_w[6])
for i in range(7, len(files_w)):
    w_section0 = np.load(files_w[i])
    w_section = np.append(w_section, w_section0, axis = 0)
#%%
PAV_b_OP = b_section[75*0:75*(0+1), :]
PAV_u_OP = u_section[75*0:75*(0+1), :]
PAV_w_OP = w_section[75*0:75*(0+1), :]
for j in range(1, 8):
    print(j)
    PAV_b_OP = PAV_b_OP + b_section[75*j:75*(j+1), :]
    PAV_u_OP = PAV_u_OP + u_section[75*j:75*(j+1), :]
    PAV_w_OP = PAV_w_OP + w_section[75*j:75*(j+1), :]
PAV_b_OP = PAV_b_OP/8
PAV_u_OP = PAV_u_OP/8
PAV_w_OP = PAV_w_OP/8
# %%
dir_out = file_dir + "/TX/"
if not os.path.exists(dir_out):
    os.makedirs(dir_out)
    
#%%
np.save(dir_out + "PAV_b_OP.npy", PAV_b_OP)
np.save(dir_out + "PAV_u_OP.npy", PAV_u_OP)
np.save(dir_out + "PAV_v_OP.npy", PAV_w_OP)
# %%
