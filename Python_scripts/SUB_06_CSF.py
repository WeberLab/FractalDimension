#!/usr/bin/env python3

import multiprocessing
from joblib import Parallel, delayed
import os 
import nibabel as nib
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import FracTool_Current_Adjustments



csf_matter = 'Masks_SUB_06/csf.nii.gz'
csf_img = nib.load(csf_matter, mmap=False)
csf_array = csf_img.get_fdata()
csf_sq = csf_array[:,:,10,0:449]
#imgplot = plt.imshow(slice_sq[:,:,1])
#imgplot = plt.imshow(csf_sq[:,:,1],cmap='Greys',  interpolation='nearest')
#Grey_fig = imgplot.get_figure()
#plt.show()
[N1,N2,N3] = csf_sq.shape
row = np.arange(0,N1)
column = np.arange(0,N2)
TR = 2.0 #temporal resolution of signal 

def FracTool_voxel(i, j):
    '''This function returns the Hurst coefficient and class of each voxel in the fMRI slice'''
    global TR
    rawbold = (csf_sq[i,j])
    result = FracTool_Current_Adjustments.FracTool(rawbold, TR) #run Fractool on each signal in signal array
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

Hurst_matrix = Hurst_matrix.reshape(N1,N2)
np.savetxt('Masks_SUB_06/Results/Hurst_matrix_SUB_06_CSF_450.txt',Hurst_matrix,fmt='%.2f')
Hurst_map = sns.heatmap(Hurst_matrix)
Hurst_fig = Hurst_map.get_figure()
Hurst_fig.suptitle('Hurst Coefficient')
Hurst_fig.savefig('Masks_SUB_06/Results/Hurst_Heatmap_SUB_06_CSF_450.png')
plt.close(Hurst_fig)


Class_matrix = Class_matrix.reshape(N1,N2)
np.savetxt('Masks_SUB_06/Results/Class_matrix_SUB_06_CSF_450.txt',Class_matrix,fmt='%.2f')
Class_map = sns.heatmap(Class_matrix, vmin=0, vmax=3)
Class_fig = Class_map.get_figure()
Class_fig.suptitle('Class')
Class_fig.savefig('Masks_SUB_06/Results/Class_Heatmap_SUB_06_CSF_450.png')
plt.close(Class_fig)