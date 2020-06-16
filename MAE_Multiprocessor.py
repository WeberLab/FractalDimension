#!/usr/bin/env python3

import multiprocessing
import numpy as np
from joblib import Parallel, delayed
#from ipynb.fs.defs.Fractool_Current import FracTool
import FracTool_Current
import pickle

#open saved expected Hurst array and signal array
with open("Hurst_expected_pickle.txt", "rb") as fp:
    Hurst_expected = pickle.load(fp)
with open("signal_pickle.txt", "rb") as fp:
    raw_sig = pickle.load(fp)
raw_sig = np.asarray(raw_sig) #convert to numpy array
n = len(Hurst_expected) #number of signals
output = range(n) #initialize results array
#generate expected class array
Class_expected = np.tile([2,0], n) #since signals are alternating fGm/fBm, class will alternate between 1 and 2
#initialize empty arrays for actual Hurst and class results from FracTool
Hurst_actual = np.zeros(n) 
Class_actual = np.zeros(n)
misclassified_total = 0
diff = 0
summation = 0
def MAE_calc(i):
    '''This function returns the MAE and miscalculation rate (MR) of FracTool, with comparision to known Hurst values 
    of the generated signal array of both fbm and fgn signals'''
    global misclassified_total #allows variables to exist in function 
    global diff
    result = FracTool_Current.FracTool(raw_sig[i]) #run Fractool on each signal in signal array
    Hurst_actual[i] = result[1] #result[1] of FracTool is Hurst value
    Class_actual[i] = result[0] #result[0] of FracTool is Class
    if Class_actual[i] != Class_expected[i]: #if signal is misclassified
        misclassified_total = 1 #increase misclassified total by 1
        Hurst_expected[i] = 0 #set Hurst expected and actual to 0 - excludes them from MAE calculation
        Hurst_actual[i] = 0
    diff = abs(Hurst_actual[i] - Hurst_expected[i])
    return [misclassified_total, diff]

#multiprocessing for loop
num_cores = multiprocessing.cpu_count()
output = Parallel(n_jobs=num_cores)(delayed(MAE_calc)(i) for i in output)
#output is array of the signal array size with two columns
#column one is 0 or 1 (not misclassified or misclassified) and column two is the difference between the signal's
#expected Hurst and the actual Hurst FracTool generated
output = np.array(output)
#print(output)
misclassified_total = np.sum(output[:,0]) #sum of column one is the total number of signals misclassified
summation = np.sum(output[:,1]) #sum of column 2 is the total difference of all signals
MAE = summation/n #MAE is summation divided by total number of signals
MR = misclassified_total/n #MR is total number of misclassified signals divided by number of signals
print(MAE, MR)