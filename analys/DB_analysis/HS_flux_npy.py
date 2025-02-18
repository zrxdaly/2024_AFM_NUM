# %%
import numpy as np
import glob
import os 
import sys

file_dir = os.path.dirname(os.path.realpath(__file__))
# %%
# for idx in np.arange(10):
#     heig = str(np.arange(5, 105, 10)[idx]).zfill(3)
#     filesB_HS = sorted(glob.glob(file_dir + "/HS/HS_By%s*.npy"%heig))
#     filesU_HS = sorted(glob.glob(file_dir + "/HS/HS_Uy%s*.npy"%heig))
#     filesW_HS = sorted(glob.glob(file_dir + "/HS/HS_Wy%s*.npy"%heig))
#     # % the array now is t, y, z
#     for i in range(len(filesB_HS)):
#         sqe = str(i).zfill(2)   
#         print(i)
#         B_HS = np.load(filesB_HS[i])
#         U_HS = np.load(filesU_HS[i])
#         W_HS = np.load(filesW_HS[i])
#         Bf_HS = B_HS - np.mean(B_HS, axis = (1,2))[:, np.newaxis, np.newaxis]
#         Uf_HS = U_HS - np.mean(U_HS, axis = (1,2))[:, np.newaxis, np.newaxis]
#         Wf_HS = W_HS - np.mean(W_HS, axis = (1,2))[:, np.newaxis, np.newaxis]
#         BW_flux = Bf_HS * Wf_HS
#         UW_flux = Uf_HS * Wf_HS
#         np.save(file_dir + "/HS/BW_flux_y%s_%s.npy"%(heig, sqe), BW_flux)
#         np.save(file_dir + "/HS/UW_flux_y%s_%s.npy"%(heig, sqe), UW_flux)

# %%
idx = int (sys.argv[1])
heig = str(np.arange(5, 105, 10)[idx]).zfill(3)
filesB_HS = sorted(glob.glob(file_dir + "/HS/HS_By%s*.npy"%heig))
filesU_HS = sorted(glob.glob(file_dir + "/HS/HS_Uy%s*.npy"%heig))
filesW_HS = sorted(glob.glob(file_dir + "/HS/HS_Wy%s*.npy"%heig))
# % the array now is t, y, z
for i in range(len(filesB_HS)):
    sqe = str(i).zfill(2)
    print(i)
    B_HS = np.load(filesB_HS[i])
    U_HS = np.load(filesU_HS[i])
    W_HS = np.load(filesW_HS[i])
    Bf_HS = B_HS - np.mean(B_HS, axis = (1,2))[:, np.newaxis, np.newaxis]
    Uf_HS = U_HS - np.mean(U_HS, axis = (1,2))[:, np.newaxis, np.newaxis]
    Wf_HS = W_HS - np.mean(W_HS, axis = (1,2))[:, np.newaxis, np.newaxis]
    BW_flux = Bf_HS * Wf_HS
    UW_flux = Uf_HS * Wf_HS
    np.save(file_dir + "/HS/BW_flux_y%s_%s.npy"%(heig, sqe), BW_flux)
    np.save(file_dir + "/HS/UW_flux_y%s_%s.npy"%(heig, sqe), UW_flux)