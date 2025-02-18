# %%
import numpy as np
import glob
import os 
import sys

# %%
idx = int (sys.argv[1])
heig = str(np.arange(5, 30, 10)[idx]).zfill(3)
file_dir = os.path.dirname(os.path.realpath(__file__))
files_HS = sorted(glob.glob(file_dir + "/HS/HS_TVy%s*.npy"%heig))
vaR = files_HS[0][-13:-11]
# %% the array now is t, y, z
RE_B = np.load(files_HS[0])
for i in range(1, 6):
    RE_B0 = np.load(files_HS[i])
    RE_B = RE_B + RE_B0
RE_B = RE_B/6
np.save(file_dir + "/HS/PAV_TV%s_RE.npy"%heig, RE_B)

del RE_B
del RE_B0
# %% the array now is t, y, z
OP_B = np.load(files_HS[6])
for i in range(7, len(files_HS)):
    OP_B0 = np.load(files_HS[i])
    OP_B = np.append(OP_B, OP_B0, axis = 0)
    
PAV_OP_B = OP_B[72*0:72*(0+1), :, :]
for i in range(1, 7):
    print(i)
    PAV_OP_B = PAV_OP_B + OP_B[72*i:72*(i+1), :, :]
    
PAV_OP_B = PAV_OP_B/7
np.save(file_dir + "/HS/PAV_TV%s_OP.npy"%heig, PAV_OP_B)
