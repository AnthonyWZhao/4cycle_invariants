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

try:
    opts, args = getopt.getopt(sys.argv[1:], "hi:o:s:e:")
except getopt.GetoptError:
    print("Option not recognised") 
    print("python read_plot.py -i <input file> -o <output file>")
    print("python read_plot.py -h for further guidance")
    sys.exit(2)
for opt, arg in opts:
    if opt == "-h":
        print("python read_plot.py -i <input file> -o <output file>")
        print("-i <input file>")
        print("-o <output file>")
        sys.exit()
    elif opt in ("-i"):
        inputFilename = arg
    elif opt in ("-o"):
        outputFilename = arg
    elif opt in ("-s"):
        start = int(arg)
    elif opt in ("-e"):
        end = int(arg) + 1       
    else:
        print("ERROR: Could not understand option " + opt)
        sys.exit(2)

if len(inputFilename) == 0:
    print("ERROR: You must provide an input file with -i")
    sys.exit(2)
#if len(outputFilename) == 0:
    #print("ERROR: You must provide an output file with -o")
    #sys.exit(2)

 
combinedDataDict = defaultdict(list) 

for i in range(start, end):   
    with open(inputFilename + "_" + str(i) + ".txt") as f:
        for line in f:
            temp_line = line.split()[1::]
            temp_line[1] = float(temp_line[1]) 
            #remove 'Taxon'
            temp_line[0] = temp_line[0].replace('Taxon', '') 
            combinedDataDict[temp_line[0]].append(temp_line[1])

#remove largest 
combinedDataDict.pop('(1,2,0,3)', None)
combinedDataDict.pop('(2,0,3,1)', None)
combinedDataDict.pop('(0,2,1,3)', None)
combinedDataDict.pop('(3,0,2,1)', None)

fig, ax = plt.subplots()

ax.boxplot(combinedDataDict.values(), vert=False)

ax.set_yticklabels(combinedDataDict.keys())                

plt.show() 
                
        
        