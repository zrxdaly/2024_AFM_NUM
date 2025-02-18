#%%
import numpy as np
import glob
import os 
import sys

# %%
# idx = int (sys.argv[1])
# heig = str(np.arange(5, 40, 10)[idx]).zfill(3)
file_dir = os.path.dirname(os.path.realpath(__file__))
dir_out = file_dir + "/HS/"
if not os.path.exists(dir_out):
    os.makedirs(dir_out)
    
for idx in np.arange(10):
    heig = str(np.arange(5, 105, 10)[idx]).zfill(3)
    # files_U = sorted(glob.glob(file_dir[:-6] + "u_slice/U%st=*"%heig))
    # files_V = sorted(glob.glob(file_dir[:-6] + "u_slice/V%st=*"%heig))
    files_W = sorted(glob.glob(file_dir[:-6] + "u_slice/W%st=*"%heig))
    vaR = "W"
    # the array now is t, y, z
    for i in range(0, len(files_W)):
        sqe = str(i).zfill(2)
        print(i)
        # USlice = np.fromfile(files_U[i], dtype = float).reshape(25, 2049, 2049)
        # VSlice = np.fromfile(files_V[i], dtype = float).reshape(25, 2049, 2049)
        WSlice = np.fromfile(files_W[i], dtype = float).reshape(25, 2049, 2049)
        # MSlice = (USlice**2 + VSlice**2 + WSlice**2)**0.5
        np.save(dir_out + "HS_%sy%s_%s.npy"%(vaR, heig, sqe), WSlice)
    
ls # %%
