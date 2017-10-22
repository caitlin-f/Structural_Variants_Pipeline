#!/usr/bin/env python3
"""
Caitlin Falconer - 22/10/2017

Called within largeSV_pipeline.sh.
Extract relevant columns from output of BreakDancer, Pindel, Crest and Delly
Export to file:
# Sample name:		Reference:
SVType		StartPos	Length		Tool

Usage: python3 collate_output.py Path/To/Output/2_LargeSVs SampleFilename ReferenceFilename
"""

import csv
import os
import sys
csv.field_size_limit(sys.maxsize)

DIR = sys.argv[1] 
FILE = sys.argv[2]
REF= sys.argv[3]
output = open("{}/{}.txt".format(DIR, FILE), 'a')
output.write("##Sample: {}\tReference: {}\n".format(FILE, REF))
output.write("#SVType\tStartPos\tSize\tTool\n")

# BreakDancer
try:
	with open("{}/BreakDancer/{}.bd.out".format(DIR, FILE)) as tsv:
		for line in csv.reader(tsv, delimiter="\t"):
			if line[0].startswith('#'):
				pass
			else:
				output.write("{}\t{}\t{}\tBreakDancer-max\n".format(line[6],line[1],line[7]))
except FileNotFoundError:
	pass

# Crest
# .predSV.txt output format:
# left_chr  left_pos  left_strand  num_left_soft-clipped reads 
# right_chr  right_pos  right_strand  num_right_soft-clipped reads
# SV type  ...
try:
	with open("{}/Crest/{}.bam.predSV.txt".format(DIR, FILE)) as tsv:
		for line in csv.reader(tsv, delimiter="\t"):
			output.write("{}\t{}\t{}\tCrest\n".format(line[8],line[1],int(line[5])-int(line[1])))
except FileNotFoundError:
	pass


# Delly
try:
	with open("{}/Delly/{}.vcf".format(DIR, FILE)) as tsv:
		for line in csv.reader(tsv, delimiter="\t"):
			if line[0].startswith('#'):
				pass
			else:
				for info in line[7].split(";"):
					if info.startswith('SVTYPE'):
						SV=info.split("=")[1]
					if info.startswith('END'):
						END=info.split("=")[1]
				output.write("{}\t{}\t{}\tDelly\n".format(SV,line[1],int(END)-int(line[1])))
except FileNotFoundError:
	pass


# Pindel
try:
	for filename in os.listdir("{}/Pindel/vcf".format(DIR)):
		# Do not include Pindel short insertion events (only long insertions)
		if (filename.startswith("{}".format(FILE))) and ('SI' not in filename): # Not collecting short-intersions
			with open("{}/Pindel/vcf/{}".format(DIR, filename)) as tsv:
				for line in csv.reader(tsv, delimiter="\t"):
					if line[0].startswith('#'):
						pass
					else:
						for n,info in enumerate(line[7].split(";")):
							if info.startswith('SVTYPE'):
								SV=info.split("=")[1]
							if info.startswith('SVLEN'):
								LEN=int(info.split("=")[1])
								# Swap sign
								if LEN > 0:
									LEN *= -1
								elif LEN < 0:
									LEN = abs(LEN)
						if 'q' in filename:
							output.write("{}\t{}\t{}\tPindel-q\n".format(SV,line[1],LEN))
						else:
							output.write("{}\t{}\t{}\tPindel\n".format(SV,line[1],LEN))
except FileNotFoundError:
	pass
	
# Dispersed duplication file does not convert to vcf properly
try:
	with open("{}/Pindel/{}_Results/{}.pd.q.out_DD".format(DIR, FILE, FILE)) as tsv:
		for line in csv.reader(tsv, delimiter="\t"):
			if line[0].startswith('#'):
				pass
			elif line[1] == 'DD' :
				output.write("DDup\t{}\t?\tPindel-q\n".format(line[4]))
except FileNotFoundError:
	pass


