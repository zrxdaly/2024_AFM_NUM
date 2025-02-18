# %%
import numpy as np
import glob
import os

# %%
file_dir = os.path.dirname(os.path.realpath(__file__))
files_Vali = sorted(glob.glob(file_dir[:-6] + "VALI_sim/b2_t=*"))
# %% the array now is t, y, z
Slice = np.loadtxt(files_Vali[0]).reshape(25, 311, 311)
for i in range(1, len(files_Vali)):
    Slice0 = np.loadtxt(files_Vali[i]).reshape(25, 311, 311)
    Slice = np.append(Slice, Slice0, axis = 0)
#%%
Slice = Slice/9.81*273.15
#%% check the variable
import matplotlib.pyplot as plt
plt.plot(Slice[:,100,100])
# %% save the variable
dir_out = file_dir + "/W12/"
if not os.path.exists(dir_out):
    os.makedirs(dir_out)
#%%
np.save(dir_out + "b2.npy", Slice)
# %%
