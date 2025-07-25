#!/usr/bin/env python

import sys, getopt, errno, time
import numpy as np
import mpmath
from skbio import DNA, TabularMSA
from multiprocessing import Pool
import utils.polynomials
from utils.polynomials import Polynomial
import operator
import random
from collections import defaultdict
import matplotlib.pylab as plt
import ruptures as rpt
import statistics

#Metric Functions:
def get1Norm(values):
	total = mpmath.mpf(0)
	for val in values:
		total = mpmath.fadd(total, mpmath.fabs(val))
	return total

def get1NormNormalised(values):
	total = get1Norm(values)
	return total/len(values)

def getMaxNorm(values):
	maxVal = 0
	for val in values:
		if abs(val) > maxVal:
			maxVal = abs(val)
	return maxVal

def getEuclideanNorm(values):
	total = 0
	for val in values:
		total += val*val
	return np.sqrt(total)

def getProdScore(values):
	return mpmath.fprod(values)

def getClimbingScore(networkString):
	score = 1
	for poly in networkPolys:
		polyValues = invariant_values[poly.getPolyString()]
		sortedPolyValues = sorted(polyValues.items(), key=operator.itemgetter(1))
		for i in range(0,len(sortedPolyValues)):
			if sortedPolyValues[i][0] == networkString:
				score *= i+1
				break
	return score

Nucl = {
  "A": 0, # = (0,0)
  "C": 1, # = (0,1)
  "G": 2, # = (1,0)
  "T": 3  # = (1,1)
}

# Function to return character values when calculating Fourier transform
def Chi(char, nucl):
	if char == 0: # Chi_A
		return 1
	elif char == 1: # Chi_C
		if nucl == 1 or nucl == 3:
			return -1
		else:
			return 1
	elif char == 2: # Chi_G
		if nucl == 2 or nucl == 3:
			return -1
		else:
			return 1
	else: # Chi_T
		if nucl == 1 or nucl == 2:
			return -1
		else:
			return 1

def getFrequenciesFromPermutation(frequencyArray, permutation):
	returnArray = np.zeros(shape=(4,4,4,4), dtype=np.longdouble)
	for i in range(4):
		for j in range(4):
			for k in range(4):
				for l in range(4):
					index=(i,j,k,l)
					permutedIndex=(index[permutation[0]], index[permutation[1]], index[permutation[2]], index[permutation[3]])
					returnArray[permutedIndex] = frequencyArray[index]
	return returnArray

def evaluateBootstrap(MSA, invariants, numSamples, seed, model, first, second):
	#print("Sampling " + str(x+1) + " of " + str(numberOfBootstraps))
	invariant_values = dict()
	permutations = [(0,1,2,3),(0,2,1,3),(0,1,3,2),(1,2,0,3),(1,0,2,3),(1,0,3,2),(2,1,0,3),(2,0,1,3),(2,0,3,1),(3,1,0,2),(3,0,1,2),(3,0,2,1)]
	random.seed(seed)
	allScores = dict() 

	for poly in invariants:
		invariant_values[poly.getPolyString()] = dict()

	#pertubation parameter likely to be changed in the future
	pertubationParameter = 0.1
	
	# Create frequencies array from MSA
	originalFrequencies = np.zeros(shape=(4,4,4,4), dtype=np.longdouble)
	count = 0
	while count < numSamples:
		i = random.randrange(MSA.shape[1])
		col = []
		for j in range(4):
			col.append(str(MSA[j,i]))
		
		

		if "-" not in col: 
			if col[first] == col[second]:
				pert = random.uniform(0,1) 
				if pert < pertubationParameter:
					possible = ["A", "G", "T", "C"] 
					possible.remove(col[second])
					col[first] = random.choice(possible)
   
			originalFrequencies[Nucl[col[0]], Nucl[col[1]], Nucl[col[2]], Nucl[col[3]]] += 1
			count += 1
 
	for i in range(4):
		for j in range(4):
			for k in range(4):
				for l in range(4):
					freq = originalFrequencies[i,j,k,l]
					if freq != 0:
						originalFrequencies[i,j,k,l] = mpmath.fdiv(freq, count)

	originalTransformed = np.zeros(shape=(4,4,4,4), dtype=np.longdouble)
	for i in range(4):
		for j in range(4):
			for k in range(4):
				for l in range(4):
					transformed_value = mpmath.mpf(0)
					for w in range(4):
						for x in range(4):
							for y in range(4):
								for z in range(4):
									if originalFrequencies[w,x,y,z] > 0:
										transformed_value = mpmath.fadd(transformed_value, mpmath.fprod([Chi(i,w), Chi(j,x), Chi(k,y), Chi(l,z), originalFrequencies[w,x,y,z]]))
					originalTransformed[i,j,k,l] = mpmath.fdiv(transformed_value, mpmath.power(4,4))

	if model == "JC":
		# Average over JC classes
		val = mpmath.fdiv(originalTransformed[0,0,1,1] + originalTransformed[0,0,2,2] + originalTransformed[0,0,3,3], 3)
		originalTransformed[0,0,1,1] = val
		originalTransformed[0,0,2,2] = val
		originalTransformed[0,0,3,3] = val

		val = mpmath.fdiv(originalTransformed[0,1,0,1] + originalTransformed[0,2,0,2] + originalTransformed[0,3,0,3], 3)
		originalTransformed[0,1,0,1] = val
		originalTransformed[0,2,0,2] = val
		originalTransformed[0,3,0,3] = val

		val = mpmath.fdiv(originalTransformed[0,1,1,0] + originalTransformed[0,2,2,0] + originalTransformed[0,3,3,0], 3)
		originalTransformed[0,1,1,0] = val
		originalTransformed[0,2,2,0] = val
		originalTransformed[0,3,3,0] = val

		val = mpmath.fdiv(originalTransformed[0,1,2,3] + originalTransformed[0,1,3,2] + originalTransformed[0,2,1,3] +
											originalTransformed[0,2,3,1] + originalTransformed[0,3,1,2] + originalTransformed[0,3,2,1], 6)
		originalTransformed[0,1,2,3] = val
		originalTransformed[0,1,3,2] = val
		originalTransformed[0,2,1,3] = val
		originalTransformed[0,2,3,1] = val
		originalTransformed[0,3,1,2] = val
		originalTransformed[0,3,2,1] = val

		val = mpmath.fdiv(originalTransformed[1,0,0,1] + originalTransformed[2,0,0,2] + originalTransformed[3,0,0,3], 3)
		originalTransformed[1,0,0,1] = val
		originalTransformed[2,0,0,2] = val
		originalTransformed[3,0,0,3] = val

		val = mpmath.fdiv(originalTransformed[1,0,2,3] + originalTransformed[1,0,3,2] + originalTransformed[2,0,1,3] +
											originalTransformed[2,0,3,1] + originalTransformed[3,0,1,2] + originalTransformed[3,0,2,1], 6)
		originalTransformed[1,0,2,3] = val
		originalTransformed[1,0,3,2] = val
		originalTransformed[2,0,1,3] = val
		originalTransformed[2,0,3,1] = val
		originalTransformed[3,0,1,2] = val
		originalTransformed[3,0,2,1] = val

		val = mpmath.fdiv(originalTransformed[1,1,0,0] + originalTransformed[2,2,0,0] + originalTransformed[3,3,0,0], 3)
		originalTransformed[1,1,0,0] = val
		originalTransformed[2,2,0,0] = val
		originalTransformed[3,3,0,0] = val

		val = mpmath.fdiv(originalTransformed[1,1,1,1] + originalTransformed[2,2,2,2] + originalTransformed[3,3,3,3], 3)
		originalTransformed[1,1,1,1] = val
		originalTransformed[2,2,2,2] = val
		originalTransformed[3,3,3,3] = val

		val = mpmath.fdiv(originalTransformed[1,2,0,3] + originalTransformed[1,3,0,2] + originalTransformed[2,1,0,3] +
											originalTransformed[2,3,0,1] + originalTransformed[3,1,0,2] + originalTransformed[3,2,0,1], 6)
		originalTransformed[1,2,0,3] = val
		originalTransformed[1,3,0,2] = val
		originalTransformed[2,1,0,3] = val
		originalTransformed[2,3,0,1] = val
		originalTransformed[3,1,0,2] = val
		originalTransformed[3,2,0,1] = val

		val = mpmath.fdiv(originalTransformed[1,2,1,2] + originalTransformed[1,3,1,3] + originalTransformed[2,1,2,1] +
											originalTransformed[2,3,2,3] + originalTransformed[3,1,3,1] + originalTransformed[3,2,3,2], 6)
		originalTransformed[1,2,1,2] = val
		originalTransformed[1,3,1,3] = val
		originalTransformed[2,1,2,1] = val
		originalTransformed[2,3,2,3] = val
		originalTransformed[3,1,3,1] = val
		originalTransformed[3,2,3,2] = val

		val = mpmath.fdiv(originalTransformed[1,2,3,0] + originalTransformed[1,3,2,0] + originalTransformed[2,1,3,0] +
											originalTransformed[2,3,1,0] + originalTransformed[3,1,2,0] + originalTransformed[3,2,1,0], 6)
		originalTransformed[1,2,3,0] = val
		originalTransformed[1,3,2,0] = val
		originalTransformed[2,1,3,0] = val
		originalTransformed[2,3,1,0] = val
		originalTransformed[3,1,2,0] = val
		originalTransformed[3,2,1,0] = val

		val = mpmath.fdiv(originalTransformed[1,1,2,2] + originalTransformed[1,1,3,3] + originalTransformed[2,2,1,1] +
											originalTransformed[2,2,3,3] + originalTransformed[3,3,1,1] + originalTransformed[3,3,2,2], 6)
		originalTransformed[1,1,2,2] = val
		originalTransformed[1,1,3,3] = val
		originalTransformed[2,2,1,1] = val
		originalTransformed[2,2,3,3] = val
		originalTransformed[3,3,1,1] = val
		originalTransformed[3,3,2,2] = val

		val = mpmath.fdiv(originalTransformed[1,2,2,1] + originalTransformed[1,3,3,1] + originalTransformed[2,1,1,2] +
											originalTransformed[2,3,3,2] + originalTransformed[3,1,1,3] + originalTransformed[3,2,2,3], 6)
		originalTransformed[1,2,2,1] = val
		originalTransformed[1,3,3,1] = val
		originalTransformed[2,1,1,2] = val
		originalTransformed[2,3,3,2] = val
		originalTransformed[3,1,1,3] = val
		originalTransformed[3,2,2,3] = val


	elif model == "K2P":
		val = mpmath.fdiv(originalTransformed[0,0,1,1] + originalTransformed[0,0,3,3], 2)
		originalTransformed[0,0,1,1] = val
		originalTransformed[0,0,3,3] = val

		val = mpmath.fdiv(originalTransformed[0,1,0,1] + originalTransformed[0,3,0,3], 2)
		originalTransformed[0,1,0,1] = val
		originalTransformed[0,3,0,3] = val

		val = mpmath.fdiv(originalTransformed[0,1,1,0] + originalTransformed[0,3,3,0], 2)
		originalTransformed[0,1,1,0] = val
		originalTransformed[0,3,3,0] = val

		val = mpmath.fdiv(originalTransformed[0,1,2,3] + originalTransformed[0,3,2,1], 2)
		originalTransformed[0,1,2,3] = val
		originalTransformed[0,3,2,1] = val

		val = mpmath.fdiv(originalTransformed[0,1,3,2] + originalTransformed[0,3,1,2], 2)
		originalTransformed[0,1,3,2] = val
		originalTransformed[0,3,1,2] = val
	
		val = mpmath.fdiv(originalTransformed[0,2,1,3] + originalTransformed[0,2,3,1], 2)
		originalTransformed[0,2,1,3] = val
		originalTransformed[0,2,3,1] = val

		val = mpmath.fdiv(originalTransformed[1,0,0,1] + originalTransformed[3,0,0,3], 2)
		originalTransformed[1,0,0,1] = val
		originalTransformed[3,0,0,3] = val

		val = mpmath.fdiv(originalTransformed[1,0,1,0] + originalTransformed[3,0,3,0], 2)
		originalTransformed[1,0,1,0] = val
		originalTransformed[3,0,3,0] = val

		val = mpmath.fdiv(originalTransformed[1,0,2,3] + originalTransformed[3,0,2,1], 2)
		originalTransformed[1,0,2,3] = val
		originalTransformed[3,0,2,1] = val

		val = mpmath.fdiv(originalTransformed[1,0,3,2] + originalTransformed[3,0,1,2], 2)
		originalTransformed[1,0,3,2] = val
		originalTransformed[3,0,1,2] = val

		val = mpmath.fdiv(originalTransformed[1,1,0,0] + originalTransformed[3,3,0,0], 2)
		originalTransformed[1,1,0,0] = val
		originalTransformed[3,3,0,0] = val

		val = mpmath.fdiv(originalTransformed[1,1,1,1] + originalTransformed[3,3,3,3], 2)
		originalTransformed[1,1,1,1] = val
		originalTransformed[3,3,3,3] = val

		val = mpmath.fdiv(originalTransformed[1,1,2,2] + originalTransformed[3,3,2,2], 2)
		originalTransformed[1,1,2,2] = val
		originalTransformed[3,3,2,2] = val

		val = mpmath.fdiv(originalTransformed[1,1,3,3] + originalTransformed[3,3,1,1], 2)
		originalTransformed[1,1,3,3] = val
		originalTransformed[3,3,1,1] = val

		val = mpmath.fdiv(originalTransformed[1,2,0,3] + originalTransformed[3,2,0,1], 2)
		originalTransformed[1,2,0,3] = val
		originalTransformed[3,2,0,1] = val		

		val = mpmath.fdiv(originalTransformed[1,2,1,2] + originalTransformed[3,2,3,2], 2)
		originalTransformed[1,2,1,2] = val
		originalTransformed[3,2,3,2] = val

		val = mpmath.fdiv(originalTransformed[1,2,2,1] + originalTransformed[3,2,2,3], 2)
		originalTransformed[1,2,2,1] = val
		originalTransformed[3,2,2,3] = val

		val = mpmath.fdiv(originalTransformed[1,2,3,0] + originalTransformed[3,2,1,0], 2)
		originalTransformed[1,2,3,0] = val
		originalTransformed[3,2,1,0] = val

		val = mpmath.fdiv(originalTransformed[1,3,0,2] + originalTransformed[3,1,0,2], 2)
		originalTransformed[1,3,0,2] = val
		originalTransformed[3,1,0,2] = val

		val = mpmath.fdiv(originalTransformed[1,3,1,3] + originalTransformed[3,1,3,1], 2)
		originalTransformed[1,3,1,3] = val
		originalTransformed[3,1,3,1] = val

		val = mpmath.fdiv(originalTransformed[1,3,2,0] + originalTransformed[3,1,2,0], 2)
		originalTransformed[1,3,2,0] = val
		originalTransformed[3,1,2,0] = val

		val = mpmath.fdiv(originalTransformed[2,0,1,3] + originalTransformed[2,0,3,1], 2)
		originalTransformed[2,0,1,3] = val
		originalTransformed[2,0,3,1] = val

		val = mpmath.fdiv(originalTransformed[2,1,0,3] + originalTransformed[2,3,0,1], 2)
		originalTransformed[2,1,0,3] = val
		originalTransformed[2,3,0,1] = val

		val = mpmath.fdiv(originalTransformed[2,1,1,2] + originalTransformed[2,3,3,2], 2)
		originalTransformed[2,1,1,2] = val
		originalTransformed[2,3,3,2] = val

		val = mpmath.fdiv(originalTransformed[2,1,2,1] + originalTransformed[2,3,2,3], 2)
		originalTransformed[2,1,2,1] = val
		originalTransformed[2,3,2,3] = val

		val = mpmath.fdiv(originalTransformed[2,1,3,0] + originalTransformed[2,3,1,0], 2)
		originalTransformed[2,1,3,0] = val
		originalTransformed[2,3,1,0] = val

		val = mpmath.fdiv(originalTransformed[2,2,1,1] + originalTransformed[2,2,3,3], 2)
		originalTransformed[2,2,1,1] = val
		originalTransformed[2,2,3,3] = val

	for perm in permutations:
		dictionaryString = "("+ str(MSA.index[perm[0]]) + "," + str(MSA.index[perm[1]]) + "," + str(MSA.index[perm[2]]) + "," + str(MSA.index[perm[3]]) + ")"
		
		for poly in invariants:
			invariant_values[poly.getPolyString()][dictionaryString] = 0

		# The 2-cycle (24) gives the same network, so we would like the result to be the same in both cases
		# We achieve this by calculating the frequency scores for both networks.
		symmetricPerm = (perm[0], perm[3], perm[2], perm[1])
		for permutation in [perm, symmetricPerm]:

			# Create frequencies array for permuted network
			transformed = getFrequenciesFromPermutation(originalTransformed, permutation)
			for poly in invariants:
				invariantSum = mpmath.fadd(invariant_values[poly.getPolyString()][dictionaryString], mpmath.fabs(poly.evaluate_mpm(transformed)))
				invariant_values[poly.getPolyString()][dictionaryString] = invariantSum

		# divide by the number of symmetric MSAs for each invariant, to get the mean abs.
		for poly in invariants:
			invariant_values[poly.getPolyString()][dictionaryString] = mpmath.fdiv(invariant_values[poly.getPolyString()][dictionaryString], 2)

 
	for perm in permutations:
		dictionaryString = "("+ str(MSA.index[perm[0]]) + "," + str(MSA.index[perm[1]]) + "," + str(MSA.index[perm[2]]) + "," + str(MSA.index[perm[3]]) + ")"
		networkValues = list()

		for poly in invariants:
			val = invariant_values[poly.getPolyString()][dictionaryString]
			networkValues.append(val)

		score = get1Norm(networkValues)
		allScores[dictionaryString] = score
		 
	return [(key, allScores[key]) for key in allScores.keys()] 


if __name__ == '__main__':
	start = time.time()
	MSAFilename = ""
	scoringFunction = get1Norm
	climbingScore = False

	multiplierForTrees1 = 3
	multiplierForTrees2 = 2

	numberOfBootstraps = 100
	doPlots = False
	numProcesses = 4
	model = "JC"

	try:
		opts, args = getopt.getopt(sys.argv[1:],"ha:i:m:t:f:s:o:")
	except getopt.GetoptError:
		print("Option not recognised.")
		print("python evaluate_bootstrap.py -a <MSA file> -i <invariants file> -m <model> -t <threads>")
		print("python evaluate_bootstrap.py -h for further usage instructions.")
		sys.exit(2)
	for opt, arg in opts:
		if opt == "-h":
			print("python evaluate_bootstrap.py -a <MSA file> -i <invariants file> -m <model> -t <threads>")
			print("-a <MSA file>\t\t Multiple sequence alignment file.")
			print("-i <invariants file>\t\t File containing list of polynomial invariants to use in Fourier coordinates.")
			print("-m <model>\t\t Either JC or K2P.") 
			print("-t <threads>\t\t Number of threads to use.") 
			print("-o <Output Directory>\t\t Directory for storing final results")
			print("-f <First>\t\t First choice in pertubation pair")
			print("-s <Second>\t\t Second choice in pertubation pair")
			print("Edge direction from -s to -f")
			sys.exit()
		elif opt in ("-a"):
			MSAFilename = arg
		elif opt in ("-i"):
			NetworkGBFilename = arg
		elif opt in ("-m"):
			if arg.upper() == "JC":
				model = "JC"
			elif arg.upper() == "K2P":
				model = "K2P"
			else:
				print("Error: Could not understand model" + arg)
				sys.exit(2)
		elif opt in ("-t"):
			numProcesses = int(arg)
		elif opt in ("-f"):
			first = int(arg)
		elif opt in ("-s"):
			second = int(arg)
		elif opt in ("-o"):
			outputDirectory = arg
    
	if len(MSAFilename) == 0:
		print("Error: You must provide an MSA file with -a.")
		sys.exit(2)
	
	if (first not in {0,1,2,3}) or (second not in {0,1,2,3}):
		print(f"Error: Leaves are labelled 0, 1, 2, 3; ({first}, {second}) is not a valid leaf pair")
		sys.exit(2)

	try:
		originalMSA = TabularMSA.read(MSAFilename, constructor=DNA)
	except (ValueError, TypeError, skbio.io.UnrecognizedFormatError) as e:
		print(e)
		sys.exit(2)
	if not originalMSA or len(originalMSA) != 4:
		print("Error: MSA file must be a multiple sequence alignment of exactly 4 sequences.")
		sys.exit(2)

	# Read in polynomials
	#print("Reading invariants...")
	try:
		with open(NetworkGBFilename, 'r') as infile:
			line = infile.readline()
			line = line.strip("| ")
			fields = line.split()
			networkPolys = list()
			for polynomial in fields:
				if polynomial.split() != "0":
					polyObject = Polynomial(polynomial)
					if polyObject not in networkPolys:
						networkPolys.append(polyObject)
						#print("Added poly: " + polyObject.getPolyString())
	except (OSError, IOError) as e: 
		if getattr(e, 'errno', 0) == errno.ENOENT:
			print("Could not find file " + NetworkGBFilename)
		print(e)
		sys.exit(2)

	winners = dict()
	sortedDictStrings = sorted([str(originalMSA.index[0]), str(originalMSA.index[1]), str(originalMSA.index[2]), str(originalMSA.index[3])])
	winners["((" + sortedDictStrings[0] + "," + sortedDictStrings[1] + "),(" + sortedDictStrings[2] +  "," + sortedDictStrings[3] + "))"] = 0
	winners["((" + sortedDictStrings[0] + "," + sortedDictStrings[2] + "),(" + sortedDictStrings[1] +  "," + sortedDictStrings[3] + "))"] = 0
	winners["((" + sortedDictStrings[0] + "," + sortedDictStrings[3] + "),(" + sortedDictStrings[1] +  "," + sortedDictStrings[2] + "))"] = 0
	
	#Permutations of the 4-cycle network: = S_4 / <(24)>
	permutations = [(0,1,2,3),(0,2,1,3),(0,1,3,2),(1,2,0,3),(1,0,2,3),(1,0,3,2),(2,1,0,3),(2,0,1,3),(2,0,3,1),(3,1,0,2),(3,0,1,2),(3,0,2,1)]
	 

	MSALength = originalMSA.shape[1]
	numSamples =  0 
	for col in originalMSA.iter_positions(ignore_metadata=True):
		if "-" not in col:
			numSamples += 1

	print("Performing analysis on " + str(numberOfBootstraps) + " replicates...")
	results = dict()
	pool = Pool(processes=numProcesses)
	for x in range(numberOfBootstraps):
		results[x] = pool.apply_async(evaluateBootstrap, args=(originalMSA, networkPolys, numSamples, x, model, first, second))

	pool.close()
	pool.join()

	combinedScores = defaultdict(list)
	for x in range(numberOfBootstraps):
		currResults = results[x].get()
		for i in range(len(currResults)): 
			combinedScores[currResults[i][0]].append(float(currResults[i][1]))
    	  
	means = [(statistics.mean(combinedScores[k]), k) for k in combinedScores.keys()]
	sorted_means = sorted(means)

	signal = []
	for m in sorted_means:
		signal += [ float(mpmath.ln(num))  for num in combinedScores[m[1]]]
	#print(signal)    

	signal = np.array(signal)
	model = "l1"
	algo=rpt.Window(width=40, model=model).fit(signal)
	#algo = rpt.Binseg(model=model).fit(signal)
 
	#needs parameter tuning
	my_bkps = algo.predict(pen=np.log( len(signal)))
	#my_bkps = algo.predict(2)
	rpt.show.display(signal, my_bkps, figsize=(10,6))
	print(my_bkps)
	with open(f"{outputDirectory}/breakpoints_{first}_{second}.txt", "w") as f:
		for bkp in my_bkps:		
			f.write(str(bkp) + "\n")
	f.close()
      
     
	plt.savefig(f"{outputDirectory}/pertubation_plot_{first}_{second}", dpi=500)
	end = time.time()
	print("evaluate_bootstrap.py time: " + str(end - start))


