# %%
import numpy as np
import glob
import os 

# %%
var = "W"
file_dir = os.path.dirname(os.path.realpath(__file__))
files_Vali = sorted(glob.glob(file_dir[:-6] + "VALI_sim/W12_%s*x=528"%var))
# %% the array now is t, y, z
Slice = np.fromfile(files_Vali[0], dtype = float).reshape(25, 11, 11)
for i in range(1, len(files_Vali)-2):
    Slice0 = np.fromfile(files_Vali[i], dtype = float).reshape(25, 11, 11)
    Slice = np.append(Slice, Slice0, axis = 0)
#%%
import matplotlib.pyplot as plt
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
plt.plot(Slice[:,3,5])
# ax.set_ylim(-0.1,0.1)
# %% save the variable
dir_out = file_dir + "/W12/"
if not os.path.exists(dir_out):
    os.makedirs(dir_out)
#%%
np.save(dir_out + "W12_%s.npy"%var, Slice)
# %%
