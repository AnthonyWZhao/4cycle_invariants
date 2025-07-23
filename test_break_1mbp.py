import sys, getopt, errno, time
import numpy as np
import itertools
import random
import skbio
from skbio import TabularMSA, DNA
from collections import defaultdict
from itertools import product
import mpmath

start = time.time()


msa_length = 1000000
model = "JC"
#outputFilename = ""

edges = ["a", "b", "c", "d", "e", "f", "g", "h"]
nucleotides = ["A", "C", "G", "T"]
substitutionRates = dict()
generateEdgeRates = True
generateGamma = True
biologyModel = False
rate = "medium"
networkString = "(0,1,2,3)"

try:
    opts, args = getopt.getopt(sys.argv[1:], "hs:p:g:l:o:br:n:")
except getopt.GetoptError:
	print("Option not recognised.")
	print("python simulateJC3cycle.py -l <MSA length> -s <seed> -p <edge parameters> -g <gamma parameter> -o <output filename> -n <network string>")
	print("python simulateJC3cycle.py -h for further usage instructions.")
	sys.exit(2)
for opt, arg in opts:
	if opt == "-h":
		print("python simulateJC3cycle.py -l <MSA length> -s <seed> -p <edge parameters> -g <gamma parameter> -o <output filename>")
		print("-l <MSA length>\t\t Length of MSA to output")
		print("-s <seed>\t\t Integer value for seeding random numbers")
		print("-p <edge parameters>\t\t Comma-separated list of values giving the probability of mutation along each edge.")
		print("-g <gamma parameter>\t\t Floating point value giving the probability of a site evolving along edge e in the network.")
		print("-o <output filename>\t Filename for output MSA in phylip format")
		sys.exit()
	elif opt in ("-o"):
		outputFilename = arg 
	elif opt in ("-l"):
		msa_length = int(arg)
	elif opt in ("-s"):
		random.seed(int(arg))
	elif opt in ("-p"):
		parameters = arg.split(",")
		if len(parameters) != len(edges):
			print ("Error: User must supply all edge parameters or none.")
			sys.exit(2)
   
		i = 0
		for param in parameters:
			val = float(param)
			if val < 0 or val > 1:
				print ("Error: User supplied edge parameter is not in interval [0,1]: " + param)
				sys.exit(2)
			substitutionRates[edges[i]] = val
			i += 1
		generateEdgeRates = False
	elif opt in ("-g"):
		val = float(arg)
		if val < 0 or val > 1:
			print ("Error: User supplied gamma parameter is not in interval [0,1]: " + param)
			sys.exit(2)
		gamma = val
		generateGamma = False
	elif opt in ("-b"):
		biologyModel = True	
	elif opt in ("-r"):
		rate = arg
	elif opt in ("-n"):
		networkString = arg
  
networkPermutationString = networkString.strip("()").split(",")
networkPermutation = [int(leaf) for leaf in networkPermutationString]
inversePermutation = np.argsort(networkPermutation)

if generateGamma:
	gamma = random.random()
if generateEdgeRates:
	if rate == "mixed":
		substitutionRates["a"] = random.uniform(0,0.001)
		substitutionRates["b"] = random.uniform(0.001,0.01)
		substitutionRates["c"] = random.uniform(0.01, 0.1)
		substitutionRates["d"] = random.uniform(0,0.001)
		substitutionRates["e"] = random.uniform(0,0.001)
		substitutionRates["f"] = random.uniform(0.001,0.01)
		substitutionRates["g"] = random.uniform(0.001,0.01)
		substitutionRates["h"] = random.uniform(0.01, 0.1) 
	else:
		for edge in edges:
			if rate == "very_low":
				substitutionRates[edge] = random.uniform(0,0.001)
			elif rate == "low":
				substitutionRates[edge] = random.uniform(0.001,0.01)
			elif rate == "medium":
				substitutionRates[edge] = random.uniform(0.01, 0.1)
			else:
				substitutionRates[edge] = random.random()

if biologyModel:
	for edge in edges:
		length = substitutionRates[edge]
		transitionProb = (3 - 3 * np.exp(-4 * length /3.))/4.

simType = "from transition matrices"
if biologyModel:
	simType = "from rate matrices"
print("Simulating alignments " + simType + " with parameters: ")
print("Tree ratio: " + str(gamma))
print("\nedge\talpha")
print("-------------------------------")
for edge in edges:
	edge_params = str(edge) + "\t" + str(substitutionRates[edge])
	print(edge_params)
 
sequence_0 = ""
sequence_1 = ""
sequence_2 = ""
sequence_3 = ""

def mutate_JC(start_nucl, edge_param):
	rand = random.random()
	if rand > edge_param:
		return start_nucl
	else:
		nucl = ["A", "C", "G", "T"]
		nucl.remove(start_nucl)
		return nucl[random.randrange(3)]

columnFreqDictLists = [defaultdict(int)] 
count = 0

for i in range(msa_length):
    count += 1
    if count > 10000:
        columnFreqDictLists.append(defaultdict(int))
        count = 1  
    #start with uniform distribution at the root
    rv_root = nucleotides[random.randrange(4)]
    
    rv_4 = mutate_JC(rv_root, substitutionRates["h"])
    rv_5 = mutate_JC(rv_4, substitutionRates["g"]) 
    
    rv_6 = 0 #reticulation vertex
    
    if random.random() < gamma:
        rv_6 = mutate_JC(rv_4, substitutionRates["f"])
    else:
        rv_6 = mutate_JC(rv_5, substitutionRates["e"])
    
    
    column = [mutate_JC(rv_6, substitutionRates["a"])
              , mutate_JC(rv_5, substitutionRates["b"])
              , mutate_JC(rv_root, substitutionRates["c"])
              , mutate_JC(rv_root, substitutionRates["d"])]
     
    column = [column[inversePermutation[0]]
              , column[inversePermutation[1]]
              , column[inversePermutation[2]]
              , column[inversePermutation[3]]]
    
    columnFreqDictLists[-1][ ''.join(column)] += 1
 

lppLists = [defaultdict(np.longdouble) for d in columnFreqDictLists]



leafPatternCombs = [''.join(comb) for comb in product(['A', 'G', 'T', 'C'], repeat=len(['A', 'G', 'T', 'C']))]


for i, d in enumerate(columnFreqDictLists):
    
    f = open("placeholder/breakup/"+ str(i) + "_" + outputFilename, "w")
    for leafPattern in leafPatternCombs:
        lppLists[i][leafPattern] = mpmath.fdiv(d[leafPattern], 10000) 
        f.write(leafPattern + "\t" + str(lppLists[i][leafPattern] ) + "\n")
    f.close()   

end = time.time()
print("simulate.py time: " + str(end - start))