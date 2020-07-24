#!/usr/bin/env python3

import numpy as np
import FracTool_Current
from fbm import fbm, fgn, times
import sklearn
import random as random

N1 = 6
N2 = 6
N3 = 100

Hurst_expected = np.zeros((N1,N2))
Sig_class_expected = np.zeros((N1,N2))
Signal_matrix = np.zeros((N1,N2,N3))
fgn_section = np.floor(N2/2) #end of leftmost third of columns of matrix
   
for i in range(N1):

    for j in range(0,int(fgn_section)):
    	Hurst_fgn = random.random() 
    	Hurst_expected[i][j] = Hurst_fgn
    	Sig_class_expected[i][j] = 0
    	Signal_matrix[i][j] = fgn(n=N3, hurst = Hurst_fgn, length=1, method='daviesharte')
        
    for k in range(int(fgn_section),N2):
        Hurst_fbm = random.random()
        Hurst_expected[i][k] = Hurst_fbm
        Sig_class_expected[i][k] = 2
        Signal_matrix[i][k] = fbm(n=N3+1, hurst = Hurst_fbm, length=1, method='daviesharte')[1:N3+1]
        
Hurst_actual = np.zeros((N1,N2))
Sig_class_actual = np.zeros((N1,N2))
for n in range(N1): 
    for m in range(N2):
        result = FracTool_Current.FracTool(Signal_matrix[n,m])
        Hurst_actual[n,m] = result[1] 
        Sig_class_actual[n,m] = result[0] 
    
Signal_matrix = []
Hurst_actual = np.nan_to_num(Hurst_actual)
summation = 0
misclassified_total = 0

for x in range(N1):
    for y in range(N2):
        if ((Sig_class_actual[x,y] != Sig_class_expected[x,y])):
            misclassified_total = misclassified_total + 1
            Hurst_actual[x,y] = 0
            Hurst_expected[x,y] = 0
              
                
        diff = abs((Hurst_actual[x,y] - Hurst_expected[x,y]))
        summation = summation + diff
    
Misclassified_Rate = misclassified_total/(N1*N2)
MAE = summation/(N1*N2)  

print(MAE, Misclassified_Rate)
