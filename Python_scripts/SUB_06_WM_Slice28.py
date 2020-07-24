#!/usr/bin/env python3

import multiprocessing
from joblib import Parallel, delayed
import os 
import nibabel as nib
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import FracTool_Current_2



grey_matter = 'whitematter.nii.gz'
gm_img = nib.load(grey_matter, mmap=False)
gm_array = gm_img.get_fdata()
gm_sq = gm_array[:,:,20,0:2250]
#imgplot = plt.imshow(slice_sq[:,:,1])
imgplot = plt.imshow(gm_sq[:,:,1],cmap='Greys',  interpolation='nearest')
Grey_fig = imgplot.get_figure()
plt.show()
[N1,N2,N3] = gm_sq.shape
row = np.arange(0,N1)
column = np.arange(0,N2)
TR = 2.0 #temporal resolution of signal 

def FracTool_voxel(i, j):
    '''This function returns the Hurst coefficient and class of each voxel in the fMRI slice'''
    global TR
    rawbold = (gm_sq[i,j])
    result = FracTool_Current_2.FracTool(rawbold, TR) #run Fractool on each signal in signal array
    H_val = result[1] #result[1] of FracTool is Hurst value
    Class_val = result[0] #result[0] of FracTool is Class
    return result 

#multiprocessing for loop
num_cores = multiprocessing.cpu_count()
output = Parallel(n_jobs=num_cores)(delayed(FracTool_voxel)(i,j) for i in row for j in column)

output = np.array(output)
output = output.astype(np.float64)

Class_matrix = output[:,0]
Hurst_matrix = output[:,1]

#generate heat maps

gm_mask = 'whitemattermask.nii.gz'
gm_img = nib.load(gm_mask, mmap=False)
gm_array = gm_img.get_fdata()
gm_sq = gm_array[:,:,20].flatten()

gm_mask_inv = (np.logical_not(gm_sq)).reshape(80,80)

Hurst_matrix = Hurst_matrix.reshape(N1,N2)
#np.savetxt('Masks_SUB_06/Results/Hurst_matrix_SUB_06_GM_450.txt',Hurst_matrix,fmt='%.2f')
Hurst_map = sns.heatmap(Hurst_matrix, vmin=0, vmax=1, mask = gm_mask_inv)
Hurst_fig = Hurst_map.get_figure()
Hurst_fig.suptitle('Hurst Coefficient')
Hurst_fig.savefig('SUB_06_WM_Slice20_Hurst.png')
plt.close(Hurst_fig)


Class_matrix = Class_matrix.reshape(N1,N2)
#np.savetxt('Masks_SUB_06/Results/Class_matrix_SUB_06_GM_450.txt',Class_matrix,fmt='%.2f')
Class_map = sns.heatmap(Class_matrix, vmin=0, vmax=3, mask = gm_mask_inv)
Class_fig = Class_map.get_figure()
Class_fig.suptitle('Class')
Class_fig.savefig('SUB_06_WM_Slice20_Class.png')
plt.close(Class_fig)