#!/usr/bin/env python3
#this script calculates H over a certain frequency range
#input: bold_name, freq_min, freq_max
#bold_name: the name of the .nii.gz file
#freq_min: the minimum frequency (inclusive)
#freq_max: the maximum frequency (inclusive)
#here is an example that will calculate H for 0.01-0.1Hz:
#welch_freq_sel.py bold_func.nii.gz 0.01 0.1 

import sys

if len(sys.argv) > 1:
    name = sys.argv[1]
    freq_min = sys.argv[2]
    freq_max = sys.argv[3]

else:
    name = input("Enter name:")
    freq_min = input("Enter min freq")
    freq_max = input("Enter max freq")

print(name)

import multiprocessing
from joblib import Parallel, delayed
import os 
import nibabel as nib
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
#import FracTool_Final
#import nolds 
import scipy
import math
from scipy.signal import welch as welch
from sklearn.linear_model import LinearRegression

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
freq_min = float(freq_min)
freq_max = float(freq_max)

def welch_voxel(i,j):
    global TR 
    voxel = (slice_sq[i,j])
    if np.mean(voxel) == 0:
        return None
    else:
        nperseg = math.floor(len(voxel)/8)
        noverlap = math.floor(nperseg/2)
        w = welch(voxel, fs = 1/TR, nperseg = nperseg, noverlap = noverlap)
        np.seterr(divide = 'ignore')
        freq=w[0]
        power=w[1]
        min_in = np.where(freq<freq_min)
        min_in = min_in[0][-1] + 1
        max_in = np.where(freq>freq_max)
        max_in = max_in[0][0] 
        freq=freq[min_in:max_in]
        power=power[min_in:max_in]
        x = np.log10(freq)[1:].reshape((-1, 1))
        y = np.log10(power)[1:]
        model = LinearRegression().fit(x, y)
        negbeta = model.coef_
        beta = negbeta*-1
        H = (beta + 1)/2
        return H


#multiprocessing for loop
num_cores = multiprocessing.cpu_count()
#base will be the bottom of stacking 2D slice arrays on top to create 3D image
base_welch = np.zeros((N1,N2))
#for loop iterates through every slice of brain and stacks slices on top of each other
for x in range(slice_num):

    slice_sq = slice_array[:,:,x,:]

    output_welch = Parallel(n_jobs=num_cores)(delayed(welch_voxel)(i,j) for i in row for j in column)
    output_welch = np.array(output_welch, dtype=object) #convert output into numpy array
    output_welch = output_welch.astype(np.float64) #convert type object list elements to type float
    Hurst_matrix_welch = output_welch.reshape(N1,N2)
    base_welch = np.dstack((base_welch,Hurst_matrix_welch))

name = name.rstrip("nii.gz") 
base_welch = base_welch[:,:,1:slice_num+1]
ni_img = nib.nifti1.Nifti1Pair(base_welch, None, n1_header)
nib.save(ni_img, f'{name}_Hurst_Welch.nii.gz')
