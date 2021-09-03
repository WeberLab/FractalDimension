#!/usr/bin/env python3
#outputs the PSD image (logfreq by logpower), the frequency range used, and Beta value (slope of linear regression)
#if you only want to calculate the PSD of the grey matter, use masked image as input (.nii.gz file)

import multiprocessing
from joblib import Parallel, delayed
import os 
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import scipy
import math
from scipy.signal import welch as welch
from sklearn.linear_model import LinearRegression
import sys
import nibabel as nib
import sys

if len(sys.argv) > 1:
    name = sys.argv[1]
else:
    name = input("Enter name:")

outDir=sys.argv[2]

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
row = np.arange(N1)
column = np.arange(N2)
#calculate frequency vector from single voxel 
samp_slice = slice_array[:,:,slice_num-1,:]
vox = (samp_slice[N1-1,N2-1])
nperseg = math.floor(len(vox)/8)
noverlap = math.floor(nperseg/2)
w = welch(vox, fs = 1/TR, nperseg = nperseg, noverlap = noverlap)
np.seterr(divide = 'ignore')
freq=w[0]


def PSD_voxel(i,j):
    '''calculates the power spectrum of one voxel using Welchs method'''
    global TR 
    global freq
    voxel = (slice_sq[i,j])
    if np.mean(voxel) == 0:
        return np.zeros(len(freq))
    else:
        nperseg = math.floor(len(voxel)/8)
        noverlap = math.floor(nperseg/2)
        w = welch(voxel, fs = TR, nperseg = nperseg, noverlap = noverlap)
        np.seterr(divide = 'ignore')
        power=w[1]
        return power


#multiprocessing for loop
num_cores = multiprocessing.cpu_count()

#for loop iterates through every slice of brain and stacks slices on top of each other
tot = np.zeros(len(freq))
count = 0
for x in range(slice_num):
    slice_sq = slice_array[:,:,x,:]
    out = Parallel(n_jobs=num_cores)(delayed(PSD_voxel)(i,j) for i in row for j in column)
    out = np.array(out, dtype=object) 
    for x in np.arange(len(out)):
        tot = tot + out[x]
        if np.mean(out[x]) == 0:
            pass
        else:
            count = count+1
            tot = tot + out[x]
            
#calculate the average power spectra array from all the voxels
avg=tot/count
avg = avg.astype(float)

#get x and y axis for freq and avg power
x = np.log10(freq)[1:-1].reshape((-1, 1))
y = np.log10(avg)[1:-1]


#create plot
fig, ax = plt.subplots()
fig.canvas.draw()
ax.plot(x,y)
#x axis is log10 freq (can convert to freq with np.power(10,x))
#y axis is power
model = LinearRegression().fit(x, y)

#calculate beta and plot linear regression line
negbeta = model.coef_
beta = negbeta*-1

x_new=x
y_new=model.predict(x)
ax.plot(x_new, y_new)

minfreq = np.power(10,x[0])
maxfreq = np.power(10,x[-1])

base = os.path.basename(name)
#png_name = "%a%b.png"%(outDir,base)
png_name = str(outDir) + str(base) + ".png"
plt.savefig(png_name)
print("Freq range:", minfreq[0],"-",maxfreq[0],"Hz")
print("Beta (slope):", negbeta[0])



