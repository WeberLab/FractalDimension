#!/usr/bin/env python3

import sys

if len(sys.argv) > 1:
    name = sys.argv[1]
else:
    name = input("Enter name:")

print(name)

import multiprocessing
from joblib import Parallel, delayed
import os 
import nibabel as nib
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import FracTool_Final

#load sample slice and convert to np array
example_slice = name
slice_img = nib.load(example_slice, mmap=False)
#get header of image
n1_header = slice_img.header
#get TR from header
TR = n1_header.get_zooms()
TR = np.asarray(TR)
TR = TR[3]
slice_array = slice_img.get_fdata()

[N1, N2, slice_num, N3] = slice_array.shape

row = np.arange(0,N1)
column = np.arange(0,N2)

def FracTool_voxel(i, j):
    '''This function returns the Hurst coefficient and class of each voxel in the slice'''
    global TR
    rawbold = (slice_sq[i,j])
    result = FracTool_Final.FracTool(rawbold, TR) #run Fractool on each signal in signal array
    H_val = result[1] #result[1] of FracTool is Hurst value
    #uncomment Class_val if you want a 3D image of signal class
    #Class_val = result[0] #result[0] of FracTool is Class
    return result #result is array of 2 elements

#multiprocessing for loop
num_cores = multiprocessing.cpu_count()
#base will be the bottom of stacking 2D slice arrays on top to create 3D image
base = np.zeros((N1,N2))
#for loop iterates through every slice of brain and stacks slices on top of each other
for x in range(slice_num):
	#change 4th index in slice_array to change the number of timepoints
	#ex 0:700 is first 700, : is the whole signal 
    slice_sq = slice_array[:,:,x,:]
    output = Parallel(n_jobs=num_cores)(delayed(FracTool_voxel)(i,j) for i in row for j in column)
    output = np.array(output) #convert output into numpy array
    output = output.astype(np.float64) #convert type object list elements to type float
    #Class_matrix = output[:,0] #first column of result output array of result is class
    Hurst_matrix = output[:,1] #second column of result output array of result is Hurst 
    Hurst_matrix = Hurst_matrix.reshape(N1,N2)
    #uncomment the following 3 lines if you want to save a text file of the Hurst of each slice
    #num = x
    #filename = f'Outputs/Hurst_matrix{num}.txt'
    #np.savetxt(filename,Hurst_matrix,fmt='%.2f')

    base = np.dstack((base,Hurst_matrix))

#save nifti file of base to the Outputs folder in FractalDimension   
name = name.rstrip("nii.gz") 
base = base[:,:,1:slice_num+1]
ni_img = nib.nifti1.Nifti1Pair(base, None, n1_header)
nib.save(ni_img, f'{name}_Hurst.nii.gz')


