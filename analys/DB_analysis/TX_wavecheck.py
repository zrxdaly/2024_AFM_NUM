# %%
import numpy as np
import glob
import os 

# %%
file_dir = os.path.dirname(os.path.realpath(__file__))
files_TX = sorted(glob.glob(file_dir + "/TX/TX_B_*.npy"))
vaR = files_TX[0][-8]
# %% save the variable
dir_out = file_dir + "/TX/"
if not os.path.exists(dir_out):
    os.makedirs(dir_out)
# %% the array now is t, y, z
Slice = np.load(files_TX[0])[:,:,2]
for i in range(1, len(files_TX)):
    Slice0 = np.load(files_TX[i])[:,:,2]
    Slice = np.append(Slice, Slice0, axis = 0)

np.save(dir_out + "TX_wavecheck.npy", Slice)

# %%
