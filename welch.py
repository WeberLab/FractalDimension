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
#import FracTool_Final
import nolds 
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

def FracTool_voxel(i, j):
    '''This function returns the Hurst coefficient and class of each voxel in the slice'''
    global TR
    rawbold = (slice_sq[i,j])
    result_fractool = FracTool_Final.FracTool(rawbold, TR) #run Fractool on each signal in signal array
    H_val = result_fractool[1] #result[1] of FracTool is Hurst value
    #uncomment Class_val if you want a 3D image of signal class
    #Class_val = result[0] #result[0] of FracTool is Class
    return result_fractool #result is array of 2 elements

def dfa_voxel(i,j):
    rawbold = (slice_sq[i,j])
    if np.mean(rawbold) == 0:
        return None
    else:
        result_dfa = nolds.dfa(rawbold)
        return result_dfa

def hurst_voxel(i,j):
    rawbold = (slice_sq[i,j])
    result_hurst = nolds.hurst_rs(rawbold)  
    return result_hurst

def corr_dim_voxel(i,j):
    rawbold = (slice_sq[i,j])
    result_corr_dim = nolds.corr_dim(rawbold, 1)  
    return result_corr_dim

def sampen_voxel(i,j):
    rawbold = (slice_sq[i,j])
    result_sampen = nolds.sampen(rawbold)  
    return result_sampen

def welch_voxel(i,j):
    global TR 
    voxel = (slice_sq[i,j])
    if np.mean(voxel) == 0:
        return None
    else:
        nperseg = math.floor(len(voxel)/8)
        noverlap = math.floor(nperseg/2)
        w = welch(voxel, fs = TR, nperseg = nperseg, noverlap = noverlap)
        np.seterr(divide = 'ignore')
        x = np.log10(w[0])[1:].reshape((-1, 1))
        y = np.log10(w[1])[1:]
        model = LinearRegression().fit(x, y)
        negbeta = model.coef_
        beta = negbeta*-1
        H = (beta + 1)/2
        return H


#multiprocessing for loop
num_cores = multiprocessing.cpu_count()
#base will be the bottom of stacking 2D slice arrays on top to create 3D image
base_fractool = np.zeros((N1,N2))
base_dfa = np.zeros((N1,N2))
base_hurst = np.zeros((N1,N2))
base_corr_dim = np.zeros((N1,N2))
base_sampen = np.zeros((N1,N2))

base_welch = np.zeros((N1,N2))
#for loop iterates through every slice of brain and stacks slices on top of each other
for x in range(slice_num):
	#change 4th index in slice_array to change the number of timepoints
	#ex 0:700 is first 700, : is the whole signal 
    slice_sq = slice_array[:,:,x,:]

    #output_fractool = Parallel(n_jobs=num_cores)(delayed(FracTool_voxel)(i,j) for i in row for j in column)
    #output_fractool = np.array(output_fractool) #convert output into numpy array
    #output_fractool = output_fractool.astype(np.float64) #convert type object list elements to type float
    #Hurst_matrix_fractool = output_fractool[:,1] #second column of result output array of result is Hurst 
    #Hurst_matrix_fractool = Hurst_matrix_fractool.reshape(N1,N2)
    
    #output_dfa = Parallel(n_jobs=num_cores)(delayed(dfa_voxel)(i,j) for i in row for j in column)
    #output_dfa = np.array(output_dfa) #convert output into numpy array
    #output_dfa = output_dfa.astype(np.float64) #convert type object list elements to type float
    #Hurst_matrix_dfa = output_dfa.reshape(N1,N2)

    #output_hurst = Parallel(n_jobs=num_cores)(delayed(hurst_voxel)(i,j) for i in row for j in column)
    #output_hurst = np.array(output_hurst) #convert output into numpy array
    #output_hurst = output_hurst.astype(np.float64) #convert type object list elements to type float
    #Hurst_matrix_hurst = output_hurst.reshape(N1,N2)
   
    #output_corr_dim = Parallel(n_jobs=num_cores)(delayed(corr_dim_voxel)(i,j) for i in row for j in column)
    #output_corr_dim = np.array(output_corr_dim) #convert output into numpy array
    #output_corr_dim = output_corr_dim.astype(np.float64) #convert type object list elements to type float
    #Hurst_matrix_corr_dim = output_corr_dim.reshape(N1,N2)
   
    #output_sampen = Parallel(n_jobs=num_cores)(delayed(sampen_voxel)(i,j) for i in row for j in column)
    #output_sampen = np.array(output_sampen) #convert output into numpy array
    #output_sampen = output_sampen.astype(np.float64) #convert type object list elements to type float
    #Hurst_matrix_sampen = output_sampen.reshape(N1,N2)

    output_welch = Parallel(n_jobs=num_cores)(delayed(welch_voxel)(i,j) for i in row for j in column)
    output_welch = np.array(output_welch, dtype=object) #convert output into numpy array
    output_welch = output_welch.astype(np.float64) #convert type object list elements to type float
    Hurst_matrix_welch = output_welch.reshape(N1,N2)
   

    #base_fractool = np.dstack((base_fractool,Hurst_matrix_fractool))
    #base_dfa = np.dstack((base_dfa,Hurst_matrix_dfa))
    #base_hurst = np.dstack((base_hurst,Hurst_matrix_hurst))
    #base_corr_dim = np.dstack((base_corr_dim,Hurst_matrix_corr_dim))
    #base_sampen = np.dstack((base_sampen,Hurst_matrix_sampen))
    base_welch = np.dstack((base_welch,Hurst_matrix_welch))



#save nifti file of base to the Outputs folder in FractalDimension   
#name = name.rstrip("nii.gz") 
#base_fractool = base_fractool[:,:,1:slice_num+1]
#ni_img = nib.nifti1.Nifti1Pair(base_fractool, None, n1_header)
#nib.save(ni_img, f'{name}_Hurst_Fractool.nii.gz')

#name = name.rstrip("nii.gz") 
#base_dfa = base_dfa[:,:,1:slice_num+1]
#ni_img = nib.nifti1.Nifti1Pair(base_dfa, None, n1_header)
#nib.save(ni_img, f'{name}_Hurst_DFA.nii.gz')

#name = name.rstrip("nii.gz") 
#base_hurst = base_hurst[:,:,1:slice_num+1]
#ni_img = nib.nifti1.Nifti1Pair(base_hurst, None, n1_header)
#nib.save(ni_img, f'{name}_Hurst_HurstRS.nii.gz')

#name = name.rstrip("nii.gz") 
#base_corr_dim = base_corr_dim[:,:,1:slice_num+1]
#ni_img = nib.nifti1.Nifti1Pair(base_corr_dim, None, n1_header)
#nib.save(ni_img, f'{name}_Hurst_corrdim.nii.gz')


#name = name.rstrip("nii.gz") 
#base_sampen = base_sampen[:,:,1:slice_num+1]
#ni_img = nib.nifti1.Nifti1Pair(base_sampen, None, n1_header)
#nib.save(ni_img, f'{name}_Hurst_sampen.nii.gz')

name = name.rstrip("nii.gz") 
base_welch = base_welch[:,:,1:slice_num+1]
ni_img = nib.nifti1.Nifti1Pair(base_welch, None, n1_header)
nib.save(ni_img, f'{name}_Hurst_Welch.nii.gz')
