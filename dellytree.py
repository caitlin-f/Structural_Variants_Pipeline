#!/usr/bin/env python3
"""
Caitlin Falconer - 22/10/2017
Clean up collated output to remove duplication and compute distance matrix
"""

import sys
import csv
import numpy as np
import pandas as pd

# store all SV results in dictionary results{{SVTYPE:{Pos:['sample']}}}
# 0 kept as place holder to allow for first iteration of dict.keys()
results = {'DEL':{0:[[],[]]},'DUP':{0:[[],[]]}, 'INS':{0:[[],[]]}}

samplenames=[]

DIR = sys.argv[1]

def AppendResults(SVTYPE, Pos, Size):
	""" For each SVtype, add pos and sample name to results dictionary 
	Parameters:
		SVTYPE(str): Type of structural variation found
		Pos (int): Position of event (int(line[1])) 
		Size (int): Size of event (int(line[2]))"""
	for key in results[SVTYPE].keys():
		if key == 0: # first entry
			results[SVTYPE].pop(0)
			results[SVTYPE][Pos]=[[sample],[Size]]
			return
		elif (key > Pos-50) and (key < Pos+50): # append to event
			results[SVTYPE][key][0].append(sample)
			results[SVTYPE][key][1].append(Size)
			return
	# novel deletion
	results[SVTYPE][Pos]=[[sample],[Size]]

with open("{}/all_data.txt".format(DIR)) as tsv:
	for line in csv.reader(tsv, delimiter="\t"):
		try:
			if line[0].startswith("##"): # get sample name
				sample = line[0].split('\t')[0].split(' ')[1]
				samplenames.append(sample)
			if line[0].startswith("#"):
				pass
			else: # line = [SvType, Pos, Size, Tool]
				if (line[0] == 'DEL') and (line[3] == "Delly"):
					AppendResults('DEL', int(line[1]), int(line[2]))
				if (line[0] == 'DUP') and (line[3] == "Delly"):
					AppendResults('DUP', int(line[1]), int(line[2]))
				if (line[0] == 'INS') and (line[3] == "Pindel"): # Store novel large insertions detected by pindel, report size as 1
					AppendResults('INS', int(line[1]), 1)
		except IndexError:
			pass

# distance matrix[row][col]
colnames=[]
for SV in results.keys():
	colnames.extend(results[SV].keys())

matrix = np.zeros((len(samplenames),len(colnames)))
distmat = pd.DataFrame(matrix,columns=colnames,index=samplenames)

colix=-1 # column index in matrix
for SV in results.keys(): # for each structural variant key
	for pos in results[SV].keys(): # for each position (i.e column) for that structural variant
		colix+=1
		for i,SAMPLE in enumerate(results[SV][pos][0]): # for each sample (i.e. row) recorded has having SV at that Pos
			rowix=samplenames.index(SAMPLE)
			distmat.iloc[rowix,colix]=(results[SV][pos][1][i])

# normalise columns to 1
for i,col in enumerate(colnames):
	total = sum(distmat.iloc[:,i])
	distmat.iloc[:,i]=(distmat.iloc[:,i]/total)

distmat.to_csv("{}/distmat.csv".format(DIR), sep='\t')
