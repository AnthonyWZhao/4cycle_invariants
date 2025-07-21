import sys, getopt, errno, time
import numpy as np
import mpmath
from skbio import DNA, TabularMSA
import utils.polynomials
from utils.polynomials import Polynomial
import operator
import pandas as pd
from collections import defaultdict
import matplotlib.pylab as plt
import ruptures as rpt
import statistics

 

 
combinedDataDict = defaultdict(list) 

smallest = ['(1,2,4,3)','(1,2,3,4)','(2,1,4,3)','(2,1,3,4)']
 
for i in range(2001, 2100):   
    with open("placeholder/1mbp_2_sim_score/trial_score_" + str(i) + ".txt") as f:
        for line in f:
            temp_line = line.split()[1::]
            temp_line[1] = float(temp_line[1]) 
            #remove 'Taxon'
            temp_line[0] = temp_line[0].replace('Taxon', '') 
            combinedDataDict[temp_line[0]].append(temp_line[1])
            
            
 
#remove largest 
#combinedDataDict.pop('(2,3,1,4)', None)
#combinedDataDict.pop('(1,3,2,4)', None)
#combinedDataDict.pop('(3,1,4,2)', None)
#combinedDataDict.pop('(4,1,3,2)', None) 

pc = []  
for j in range(len(combinedDataDict['(1,2,4,3)'])):
    big = []
    small = []
    total = 0
    need = 0
    for network in combinedDataDict.keys():
        curr = mpmath.ln(combinedDataDict[network][j])
        total += curr
        if network in smallest:
            need += curr
            small.append(curr)
        else:
            big.append(curr)
    #print(mpmath. )            
    pc.append(mpmath.nstr(mpmath.fdiv(need,total)))
    
print(pc)            


means = [(statistics.mean(combinedDataDict[k]), k) for k in combinedDataDict.keys()]
sorted_means = sorted(means)

signal = []
for m in sorted_means:
    signal += [ float(mpmath.ln(num))  for num in combinedDataDict[m[1]]]
#print(signal)    

signal = np.array(signal)
model = "l1"
algo=rpt.Window(width=40, model=model).fit(signal)
#algo = rpt.Binseg(model=model).fit(signal)
#my_bkps = algo.predict(pen=np.log( len(signal)))
my_bkps = algo.predict(2)
rpt.show.display(signal, my_bkps, figsize=(10,6))
print(my_bkps)
plt.show()
          

