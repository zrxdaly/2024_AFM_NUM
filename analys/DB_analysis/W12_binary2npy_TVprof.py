# %%
import numpy as np
import glob
import os 

# %%
file_dir = os.path.dirname(os.path.realpath(__file__))
files_Vali = sorted(glob.glob(file_dir[:-6] + "VALI_sim/W12_TVt=*x=528"))
# %% the array now is t, y, z
Slice = np.fromfile(files_Vali[0], dtype = float).reshape(25, 11, 11)
for i in range(1, len(files_Vali)):
    Slice0 = np.fromfile(files_Vali[i], dtype = float).reshape(25, 11, 11)
    Slice = np.append(Slice, Slice0, axis = 0)
#%% chop out the duplicate part
Slice = Slice -273.15
#%%
import matplotlib.pyplot as plt
plt.plot(Slice[:,2,2])
# %% save the variable
dir_out = file_dir + "/W12/"
if not os.path.exists(dir_out):
    os.makedirs(dir_out)
# %%
np.save(dir_out + "W12_TV.npy", Slice)
# %%
