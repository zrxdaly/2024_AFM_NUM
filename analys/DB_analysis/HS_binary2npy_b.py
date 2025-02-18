#%%
import numpy as np
import glob
import os
import sys


file_dir = os.path.dirname(os.path.realpath(__file__))
dir_out = file_dir + "/HS/"
if not os.path.exists(dir_out):
    os.makedirs(dir_out)
# %%
# idx = int (sys.argv[1])
# heig = str(np.arange(5, 95, 10)[idx]).zfill(3)
for idx in np.arange(10):
    heig = str(np.arange(5, 95, 10)[idx]).zfill(3)
    files_Sl = sorted(glob.glob(file_dir[:-6] + "b_slice/B%st=*"%heig))
    vaR = files_Sl[0][-8]
    #% the array now is t, y, z
    for i in range(0, len(files_Sl)):
        sqe = str(i).zfill(2)   
        print(i)
        Slice = np.fromfile(files_Sl[i], dtype = float).reshape(25, 2049, 2049)
        Slice = Slice/9.81*273.15
        np.save(dir_out + "HS_%sy%s_%s.npy"%(vaR, heig, sqe), Slice)
# %%
