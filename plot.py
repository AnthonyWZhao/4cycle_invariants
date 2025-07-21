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
import matplotlib.pylab as labplt
import ruptures as rpt
import statistics

""" try:
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
    #sys.exit(2) """

plt.figure(figsize=(10, 20)) 

left = 10
 
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

    ax = plt.subplot(5,2, j )
    ax.boxplot(combinedDataDict.values(), vert=False)
    ax.set_title(f"1mbp_{left}_" + str(j))
    ax.set_yticklabels(combinedDataDict.keys()) 
    
plt.savefig(f"result_data/plots/box_plots/1mbp_{left}")
#plt.show()    

 

smallest = ['(1,2,4,3)','(1,2,3,4)','(2,1,4,3)','(2,1,3,4)']
 
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
    pc = []
    for ele in range(len(combinedDataDict['(1,2,4,3)'])):
        total = 0
        need = 0
        for network in combinedDataDict.keys():
            curr = mpmath.ln(combinedDataDict[network][ele])
            total += curr
            if network in smallest:
                need += curr              
        pc.append(float(mpmath.fdiv(need,total)) )
                  
    with open(f"result_data/text_data/1mbp_{left}.txt", "w+") as a:
        a.write(f"1mbp_{left}_{j}:" + "\t" 
                + f"Pc Mean: {statistics.mean(pc)}" 
                + "\t" + f"Pc Variance: {statistics.variance(pc)}" + "\n")
   
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
    my_bkps = algo.predict(pen=np.log( len(signal)))
    #my_bkps = algo.predict(2) 
    rpt.show.display(signal, my_bkps, figsize=(10,6))  
    labplt.savefig(f"result_data/plots/rupture_plots/1mbp_{left}/1mbp_{left}_{j}"  )  
    
    
                    
            
            