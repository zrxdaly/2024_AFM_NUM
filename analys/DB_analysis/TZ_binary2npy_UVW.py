# %%
import numpy as np
import glob
import os 

# %%
file_dir = os.path.dirname(os.path.realpath(__file__))
files_TX = sorted(glob.glob(file_dir[:-6] + "u_slice/Z_Wt=*"))
vaR = files_TX[0][-5]
# %% save the variable
dir_out = file_dir + "/TZ/"
if not os.path.exists(dir_out):
    os.makedirs(dir_out)
# %% the array now is t, y, z
for i in range(0, len(files_TX)):
    print(i)
    sqe = str(i).zfill(2)
    Slice = np.fromfile(files_TX[i], dtype = float).reshape(25, 100, 2049)
    np.save(dir_out + "TZ_%s_%s.npy"%(vaR, sqe), Slice)
# %%
