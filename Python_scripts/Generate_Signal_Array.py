#!/usr/bin/env python3

import pickle
import numpy as np
from fbm import fbm, fgn, times
import random
import multiprocessing
from joblib import Parallel, delayed


n = 1000 #n is number of signals  
length = 5000 #length is the length of each signal
#initialize arrays to iterate through
signals = range(n)
Hurst_expected = np.zeros(n)
#generate array of expected Hurst values 
for x in range(n):
     Hurst_expected[x] = random.random() #generate values between 0 and 1
        
def generate_signal(j):
    '''This function generates an array of signals, alternating between fBm and fGn signals'''
    Hurst_num = Hurst_expected[j] #Hurst value used is indexed from known Hurst array
    if j % 2 == 0: #every even index is a fBm signal
        sig = fbm(n=length+1, hurst = Hurst_num, length=1, method='daviesharte')[1:length+1]
    if j % 2 != 0: #every odd index is a fGn signal 
        sig = fgn(n=length, hurst = Hurst_num, length=1, method='daviesharte')    
    return sig

#multiprocessor for loop 
num_cores = multiprocessing.cpu_count()
signals = Parallel(n_jobs=num_cores)(delayed(generate_signal)(j)for j in signals)

#save/pickle signal array 
with open("signal_pickle.txt", "wb") as fp:
    pickle.dump(signals, fp)
#save/pickle Hurst expected array
with open('Hurst_expected_pickle.txt', "wb") as fp:
    pickle.dump(Hurst_expected,fp)
#save signal array as txt file for future reference    
with open('signal_copy.txt', 'w') as f:
    for signal in signals:
        f.write("%s\n" % signal)
#save Hurst expected array as txt file for future reference        
with open('Hurst_expected_copy.txt', 'w') as f:
    for Hurst in Hurst_expected:
        f.write("%s\n" % Hurst)
