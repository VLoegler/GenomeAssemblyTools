#!/usr/bin/env python
# -*- coding: utf8 -*-
#----------------------------------------------------------------------------
# Created By  : vloegler
# Created Date: 2022/09/15
# version ='1.1'
# ---------------------------------------------------------------------------
'''
This script rename contigs according to a reference genome, if at least X kb 
of a contig aligns on a reference chromosome. If a contig aligns on several
reference chromosomes, all chromosomes name will be indicated, by decreasing 
order of alignment size. 
An additionnal argument add the alignment sizes in the contig name. 
This script uses nucmer and show-coords MUMmer4's functions. 

It's better to use a reference with transposable elements masked, for example
with RepeatMasker:
	RepeatMasker reference.fasta -lib TransposableElements.fasta -s -e ncbi
and with telomeres masked.

It takes as input :
	-r --ref: the reference genome assembly (multi fasta)
	-d --draft: draft genome assemblies (multi fasta)
	-o --output: name of output fasta file
	-k --alignmentSizeKb: Minimum alignment size to have name correspondance (in kb)
	-i --alignmentSizeInfo: Add alignment size information in contig names
	-m --mummerPath: Path to mummer function, if not in path
	-t --threads: Number of threads for nucmer

Output will be fasta with renamed contigs. 
'''
# ---------------------------------------------------------------------------
import os
import sys
import argparse
import time
import csv
from datetime import datetime
from random import randint
from Tools import *
import re
# ---------------------------------------------------------------------------
# Definitions

def rank_simple(vector):
	return sorted(range(len(vector)), key=vector.__getitem__)

def decreasing_rank_simple(vector):
	return sorted(range(len(vector)), key=vector.__getitem__)[::-1]

def getShowCoords(refPath, draftPath, prefix, mummerPath, threads):
	# Run nucmer
	#nucmerShellCommand=mummerPath+"nucmer -t "+str(threads)+" --maxmatch --prefix "+prefix+" "+refPath+" "+draftPath
	nucmerShellCommand=mummerPath+"nucmer -t "+str(threads)+" --mum --prefix "+prefix+" "+refPath+" "+draftPath
	os.system(nucmerShellCommand)
	# Run show coords
	showcoordsShellCommand=mummerPath+"show-coords -TH "+prefix+".delta > "+prefix+".coords"
	os.system(showcoordsShellCommand)
	# Remove nucmer delta file
	os.remove(prefix+".delta")

# ---------------------------------------------------------------------------

# =============
# Get arguments
# =============

# Initiate the parser
parser = argparse.ArgumentParser()
parser.add_argument("-r", "--ref", help="reference genome assembly (multi fasta)", required=True)
parser.add_argument("-d", "--draft", help="draft genome assemblies (multi fasta)", required=True)
parser.add_argument("-o", "--output", help="Name of the outputed VCF", required=True)
parser.add_argument("-k", "--alignmentSizeKb", help="Minimum alignment size to have name correspondance (in kb, default 30)", type=int, default=30)
parser.add_argument("-i", "--alignmentSizeInfo", help="Add alignment size information in contig names", action='store_true')
parser.add_argument("-m", "--mummerPath", help="Path to mummer function, if not in path", type=str, default="")
parser.add_argument("-t", "--threads", help="Number of threads for nucmer", type=int, default=1)


# Read arguments from the command line
args = parser.parse_args()

refPath=args.ref
draftPath=args.draft
outputPath=args.output
alignmentSizeKb=args.alignmentSizeKb
mummer=args.mummerPath
if mummer != "" and not mummer.endswith("/") :
	mummer += "/"
threads=args.threads

print("\n\t--- RENAMING CONTIGS ---\n")
print("Arguments detected:")
print("\t--ref:\t"+refPath)
print("\t--draft:\t"+draftPath)
print("\t--output:\t"+outputPath)
print("\t--alignmentSizeKb:\t"+str(alignmentSizeKb))
if args.alignmentSizeInfo:
	print("\t--alignmentSizeInfo")
if mummer != "":
	print("\t--mummerPath:\t"+mummer)
if threads != 1:
	print("\t--threads:\t"+str(threads))
print('')


# =============================
# Get reference chromosome name
refChr=[]
ref=open(refPath, 'r')
for line in ref.readlines():
	if line.startswith(">"):
		refChr += [line.strip().split(">")[1].split(" ")[0].split("\t")[0]]
ref.close()

# =================================
# Get draft chromosome name and seq
draftChr=[]
draftSeq=[]
seq=""
draft=open(draftPath, 'r')
for line in draft.readlines():
	if line.startswith(">"):
		draftChr += [line.strip().split(">")[1].split(" ")[0].split("\t")[0]]
		if seq != "":
			draftSeq += [seq]
		seq=""
	else:
		seq += line.strip()
draftSeq += [seq]
draft.close()
draftLen=[len(x) for x in draftSeq]

refName=refPath.split("/")[-1]
draftName=draftPath.split("/")[-1]
prefix = "Alignment_" + refName + "_" + draftName + "_" + datetime.now().strftime("%Y_%m_%d_%H_%M_%S_%m_%f") + "_" + str(randint(0, 10000))

# ===============================================
# Align and get alignment coordinates with nucmer
getShowCoords(refPath, draftPath, prefix, mummer, threads)

# ==========================
# Read alignment coordinates
# Cummulative size of alignments will be stored in a matrix with 
# first dimension: Draft contig
# second dimension: Reference chromosome
x = len(draftChr)
y = len(refChr)
analysisMatrix = []
for i in range(x):
	analysisMatrix.append([])
	for j in range(y):
		analysisMatrix[i].append([])

# Fill matrix
coordsFile = open(prefix+".coords")
coords = csv.reader(coordsFile, delimiter="\t")
for row in coords:
	draftContigIndex=draftChr.index(row[8])
	refChrIndex=refChr.index(row[7])
	# Add alignment length
	chr=row[7]
	startPos=int(row[0])
	endPos=int(row[1]) + 1
	if startPos < endPos:
		coord2add = BEDcoordinates(id = chr, start = startPos, end = endPos)
	else:
		coord2add = BEDcoordinates(id = chr, start = endPos, end = endPos)
	analysisMatrix[draftContigIndex][refChrIndex] += [coord2add]
coordsFile.close()
os.remove(prefix+".coords")

# Convert all point of the matrix to BED object and get length
for i in range(len(draftChr)):
	for j in range(len(refChr)):
		analysisMatrix[i][j] = BED(analysisMatrix[i][j])

# Second matrix just containing length
x = len(draftChr)
y = len(refChr)
analysisMatrixLen = []
for i in range(x):
	analysisMatrixLen.append([])
	for j in range(y):
		analysisMatrixLen[i].append(analysisMatrix[i][j].len)


# Assign a main Chr to each contig
mainChrs = [] # Contain the main chromosome corresponding to each contig
mainChrsInfo = []
for i in range(len(draftChr)):
	mainChrs += [refChr[analysisMatrixLen[i].index(max(analysisMatrixLen[i]))]]
	mainChrsInfo += [mainChrs[i] + "=" + str(analysisMatrixLen[i][refChr.index(mainChrs[i])])]

# Check if several contigs are assigned to a single chromosome
mainChrs2 = mainChrs.copy() # Contain the main chr for each contig + indication of order if several contigs per chr
for j in range(len(refChr)):
	if mainChrs2.count(refChr[j]) > 1:
		# Get index of chromosomes
		indices = [i for i, x in enumerate(mainChrs2) if x == refChr[j]]
		# Get centers of alignments
		centers = [analysisMatrix[i][j].getCenter() for i in indices]
		suffix = [i+1 for i in rank_simple(centers)]
		for i in range(len(indices)):
			mainChrs2[indices[i]] = mainChrs2[indices[i]] + "." + str(suffix[i])

# Add secondary chromosomes (in case of translocation). That is Chromosome aligning on more than X kb on the contig
mainChrs3 = mainChrs2.copy() # Contains additional chromosome aligned on contig because of translocations
mainChrs3Info = [""] * len(draftChr)
for i in range(len(draftChr)):
	main = mainChrs[i]
	mainAlignment = analysisMatrixLen[i][refChr.index(main)]

	secondaryChrs = [x for x in refChr if x != main and analysisMatrixLen[i][refChr.index(x)] >= alignmentSizeKb*1000]
	secondaryAlignments = [analysisMatrixLen[i][refChr.index(j)] for j in secondaryChrs]
	secondaryAlignmentsOrder = decreasing_rank_simple(secondaryAlignments)

	for j in secondaryAlignmentsOrder:
		mainChrs3[i] += "_"+secondaryChrs[j]
		mainChrs3Info[i] += " "+secondaryChrs[j]+"="+str(analysisMatrixLen[i][refChr.index(secondaryChrs[j])])

# Final names
draftNewChr = []
for i in range(len(draftChr)):
	if args.alignmentSizeInfo:
		draftNewChr += [">" + mainChrs3[i] + " length=" + str(draftLen[i]) + " " + mainChrsInfo[i] + mainChrs3Info[i]]
	else:
		draftNewChr += [">" + mainChrs3[i] + " length=" + str(draftLen[i])]



# create new Fasta file with renamed contigs
out = open(outputPath, 'w')
for i in range(len(draftChr)):
	out.write(draftNewChr[i]+"\n")
	seq=re.sub("(.{80})", "\\1\n", draftSeq[i], 0, re.DOTALL)
	out.write(seq+"\n")
out.close()

print("New fasta written to: "+outputPath)
print("")

