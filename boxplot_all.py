import sys, getopt, errno, time
import numpy as np
import mpmath
from skbio import DNA, TabularMSA
import utils.polynomials
from utils.polynomials import Polynomial
import operator
import pandas as pd
from collections import defaultdict
import matplotlib.pyplot as plt


plt.figure(figsize=(20,20) ) 
for left in range(1, 11): 
    for j in range(1, 11):
        start, end = 1, 101
        combinedDataDict = defaultdict(list)
        for i in range(start, end):   
            with open(f"result_data/scores/1mbp_{left}_" + str(j) + "/score_" + str(left*1000 + (j-1)*100 + i) + ".txt") as f:
                for line in f:
                    temp_line = line.split()[1::]
                    temp_line[1] = float(temp_line[1]) 
                    #remove 'Taxon'
                    temp_line[0] = temp_line[0].replace('Taxon', '') 
                    combinedDataDict[temp_line[0]].append(temp_line[1])

        #remove largest 
        combinedDataDict.pop('(2,3,1,4)', None)
        combinedDataDict.pop('(1,3,2,4)', None)
        combinedDataDict.pop('(3,1,4,2)', None)
        combinedDataDict.pop('(4,1,3,2)', None)

        ax = plt.subplot(10,10, (left - 1)*10 + j  )
        ax.boxplot(combinedDataDict.values(), vert=False)
        ax.set_title(f"1mbp_{left}_" + str(j))
        ax.set_yticklabels(combinedDataDict.keys()) 
plt.tight_layout()    
plt.savefig(f"result_data/plots/box_plots/1mbp_all", dpi=500)
plt.show()    
