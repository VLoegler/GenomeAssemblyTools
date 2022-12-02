#!/usr/bin/env python
# -*- coding: utf8 -*-
#----------------------------------------------------------------------------
# Created By  : vloegler
# Created Date: 2022/09/28
# version ='1.3'
# ---------------------------------------------------------------------------
'''
This script rename contigs according to a reference genome, if at least X kb 
or X% of a contig aligns on a reference chromosome. If a contig aligns on 
several reference chromosomes, all chromosomes name will be indicated, by 
order of appearence in the contig. 
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
	-p --alignmentSizePercent: Minimum alignment size to have name correspondance (in % of contig size)
	-P --prefix: Prefix to add before the chromosome name
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
parser.add_argument("-k", "--alignmentSizeKb", help="Minimum alignment size to have name correspondance (in kb, default 50)", type=int, default=50)
parser.add_argument("-p", "--alignmentSizePercent", help="Minimum alignment size to have name correspondance (in percent, default 20%)", type=int, default=20)
parser.add_argument("-P", "--prefix", help="Prefix to add before the chromosome name", type=str, default="")
parser.add_argument("-i", "--alignmentSizeInfo", help="Add alignment size information in contig names", action='store_true')
parser.add_argument("-m", "--mummerPath", help="Path to mummer function, if not in path", type=str, default="")
parser.add_argument("-t", "--threads", help="Number of threads for nucmer", type=int, default=20)


# Read arguments from the command line
args = parser.parse_args()

refPath=args.ref
draftPath=args.draft
outputPath=args.output
alignmentSizeKb=args.alignmentSizeKb
alignmentSizePercent=args.alignmentSizePercent
ChrPrefix=args.prefix
mummer=args.mummerPath
if mummer != "" and not mummer.endswith("/") :
	mummer += "/"
threads=args.threads

print("\n\t--- RENAMING CONTIGS ---\n")
print("Arguments detected:")
print("\t--ref:\t"+refPath)
print("\t--draft:\t"+draftPath)
print("\t--output:\t"+outputPath)
print("\t--alignmentSizeKb:\t"+str(alignmentSizeKb)+"kb")
print("\t--alignmentSizePercent:\t"+str(alignmentSizePercent)+"%")
if ChrPrefix != "":
	print("\t--prefix:\t"+ChrPrefix)
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
refSeq=[]
seq=""
ref=open(refPath, 'r')
for line in ref.readlines():
	if line.startswith(">"):
		refChr += [line.strip().split(">")[1].split(" ")[0].split("\t")[0]]
		if seq != "":
			refSeq += [seq]
		seq=""
	else:
		seq += line.strip()
refSeq += [seq]
ref.close()
refLen=[len(x) for x in refSeq]

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
# analysisMatrixRef contains BED based on ref chr
# analysisMatrixDraft contains BED based on draft chr
x = len(draftChr)
y = len(refChr)
analysisMatrixRef = []
analysisMatrixDraft = []
for i in range(x):
	analysisMatrixRef.append([])
	analysisMatrixDraft.append([])
	for j in range(y):
		analysisMatrixRef[i].append([])
		analysisMatrixDraft[i].append([])

# Fill matrix
coordsFile = open(prefix+".coords")
coords = csv.reader(coordsFile, delimiter="\t")
for row in coords:
	draftContigIndex=draftChr.index(row[8])
	refChrIndex=refChr.index(row[7])

	# Add alignment length to analysisMatrixRef
	chr=row[7]
	startPos=int(row[0])
	endPos=int(row[1]) + 1
	if startPos < endPos:
		coord2add = BEDcoordinates(id = chr, start = startPos, end = endPos)
	else:
		coord2add = BEDcoordinates(id = chr, start = endPos, end = startPos)
	analysisMatrixRef[draftContigIndex][refChrIndex] += [coord2add]

	# Add alignment length to analysisMatrixDraft
	chr=row[8]
	startPos=int(row[2])
	endPos=int(row[3]) + 1
	if startPos < endPos:
		coord2add = BEDcoordinates(id = chr, start = startPos, end = endPos)
	else:
		coord2add = BEDcoordinates(id = chr, start = endPos, end = startPos)
	analysisMatrixDraft[draftContigIndex][refChrIndex] += [coord2add]

coordsFile.close()
os.remove(prefix+".coords")

# Convert all point of the matrix to BED object and get length
for i in range(len(draftChr)):
	for j in range(len(refChr)):
		analysisMatrixDraft[i][j] = BED(analysisMatrixDraft[i][j])
		analysisMatrixRef[i][j] = BED(analysisMatrixRef[i][j])

# Second matrix just containing length
x = len(draftChr)
y = len(refChr)
analysisMatrixLen = []
for i in range(x):
	analysisMatrixLen.append([])
	for j in range(y):
		analysisMatrixLen[i].append(analysisMatrixRef[i][j].len)

# Find correspondance between contigs and chromosomes
# A correspondance is an total alignment size greater than 30kb, 
# 20% the contig size or 20% the chromosome size
# Results will be stored in a matrix correspondanceMatrixRef[refIndex][draftIndex]
# that will contain either 0 (no correspondance between contig and chr)
# or a number which is the part of Chromosome corresponding
# 1 = first part or whole chromosome
# 2 = second part
# ...
#
# A second matrix correspondanceMatrixDraft will contain the order in which the 
# chromosome fragments are present in the contig

x = len(draftChr)
y = len(refChr)
correspondanceMatrixRef = []
correspondanceMatrixDraft = []
for i in range(x):
	correspondanceMatrixRef.append([])
	correspondanceMatrixDraft.append([])
	for j in range(y):
		correspondanceMatrixRef[i].append(0)
		correspondanceMatrixDraft[i].append(0)

# Find fragment of chromosomes
#numberFragments = [0]*len(refChr)
for j in range(len(refChr)):
	contigAssociated = [0]*len(draftChr)
	for i in range(len(draftChr)):
		sizeThreshold = min([alignmentSizeKb*1000, alignmentSizePercent*0.01*refLen[j], alignmentSizePercent*0.01*draftLen[i]])
		if analysisMatrixLen[i][j] >= sizeThreshold:
			contigAssociated[i] = analysisMatrixRef[i][j].getCenter()

	# Order chromosome fragments
	contigAssociatedOrder = [0]*len(draftChr)
	centers = [x for x in contigAssociated if x > 0]
	for i in range(len(draftChr)):
		if contigAssociated[i] != 0:
			correspondanceMatrixRef[i][j] = 1 + sorted(centers).index(contigAssociated[i])



for i in range(len(draftChr)):
	chromoAssociated = [0]*len(refChr)
	for j in range(len(refChr)):
		sizeThreshold = min([alignmentSizeKb*1000, alignmentSizePercent*0.01*refLen[j], alignmentSizePercent*0.01*draftLen[i]])
		if analysisMatrixLen[i][j] >= sizeThreshold:
			chromoAssociated[j] = analysisMatrixDraft[i][j].getCenter()

	# Order chromosome in contigs
	chromoAssociatedOrder = [0]*len(refChr)
	centers = [x for x in chromoAssociated if x > 0]
	for j in range(len(refChr)):
		if chromoAssociated[j] != 0:
			correspondanceMatrixDraft[i][j] = 1 + sorted(centers).index(chromoAssociated[j])

# Get name of each contig
draftNewChr = []
unplacedNb=0
for i in range(len(draftChr)):
	newName = ">"
	if ChrPrefix != "":
		newName += ChrPrefix+"_"

	additionnalInfo = " AlignmentSizeOnRef:"
	nbChrInContig = max(correspondanceMatrixDraft[i]) # Number of chr in the contig
	if nbChrInContig > 0:
		for n in range(nbChrInContig):
			if n != 0:
				newName += "_"
				additionnalInfo += ","
			# Add Chr Name
			j = correspondanceMatrixDraft[i].index(n+1)
			newName += refChr[j]
			# Add Chr fragment, if chromosome if fragmented
			if max([x[j] for x in correspondanceMatrixRef]) > 1:
				newName += "."+str(correspondanceMatrixRef[i][j])
			# Add alignment size
			additionnalInfo += refChr[j] + "=" + str(analysisMatrixLen[i][j])
	else:
		unplacedNb += 1
		newName += "Unplaced." + str(unplacedNb)
		additionnalInfo = ""

	if args.alignmentSizeInfo:
		draftNewChr += [newName+" length="+str(draftLen[i])+additionnalInfo]
	else:
		draftNewChr += [newName+" length="+str(draftLen[i])]


# create new Fasta file with renamed contigs
out = open(outputPath, 'w')
for i in range(len(draftChr)):
	out.write(draftNewChr[i]+"\n")
	seq=re.sub("(.{80})", "\\1\n", draftSeq[i], 0, re.DOTALL)
	out.write(seq+"\n")
out.close()

print("New fasta written to: "+outputPath)
print("")

