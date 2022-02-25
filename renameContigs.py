#!/usr/bin/env python
# -*- coding: utf8 -*-
#----------------------------------------------------------------------------
# Created By  : vloegler
# Created Date: 2022/02/24
# version ='1.0'
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

def addBEDcoordinates(b1, b2):
	if len(b1) == 0 and len(b2) == 0:
		return [[]]
	elif len(b1) == 0 :
		return [b2]
	elif len(b2) == 0:
		return [b1]
	else:
		# Get info from b1
		chr1=b1[0]
		start1=b1[1]
		end1=b1[2]
		# Get info from b2
		chr2=b2[0]
		start2=b2[1]
		end2=b2[2]
		# Store b1+b2 in newb BED object
		newb=[]
		if chr2 == chr1:
			if start2 <= start1:
				if end2 <= start1:
					newb = [b2, b1]
				elif end2 < end1:
					newb = [[chr1, start2, end1]]
				else: 
					newb = [b2] # b1 is comprised in b2
			elif start2 < end1:
				if end2 < end1:
					newb = [b1] # b2 is comprised in b1
				else: 
					newb = [[refChr, start1, end2]]
			else:
				newb = [b1, b2]
		else: 
			newb = [b1, b2]
		# Return the addition as BED object (list of BED coordinates)
		return newb

def overlapBEDcoordinates(b1, b2):
	if len(b1) == 0 and len(b2) == 0:
		return False
	elif len(b1) == 0 :
		return False
	elif len(b2) == 0:
		return False
	else:
		# Get info from b1
		chr1=b1[0]
		start1=b1[1]
		end1=b1[2]
		# Get info from b2
		chr2=b2[0]
		start2=b2[1]
		end2=b2[2]
		# Store b1+b2 in newb BED object
		newb=[]
		if chr2 == chr1:
			if start2 <= start1:
				if end2 <= start1:
					return False
				elif end2 < end1:
					return True
				else: 
					return True
			elif start2 < end1:
				if end2 < end1:
					return True
				else: 
					return True
			else:
				return False
		else: 
			return False

def addBEDoneCoordinate(BED, b):
	overlap = False
	overlappedCoordinates = []
	coordinate2pop = []
	for i in range(len(BED)):
		if overlapBEDcoordinates(BED[i], b):
			overlap = True
			overlappedCoordinates += [BED[i]]
			coordinate2pop += [i]
	if not overlap:
		return BED + [b]
	else:
		# Remove Overlapped coordinates from BED
		coordinate2pop.sort(reverse = True) # Sort reverse to remove last indexes first
		for i in coordinate2pop:
			BED.pop(i)
		# Sequentially add Overlapped coordinates to b
		for i in range(len(overlappedCoordinates)):
			b = addBEDcoordinates(b, overlappedCoordinates[i])[0]
		return BED + [b]

def addBED(BED1, BED2):
	if len(BED1) == 0 and len(BED2) == 0:
		return []
	elif len(BED1) == 0:
		return BED2
	elif len(BED2) == 0:
		return BED1
	if len(BED2) > 1:
		b = BED2.pop()
		addBED(addBED(BED1, [b]), BED2)
	else : # There is a single in BED2
		return addBEDoneCoordinate(BED1, BED2[0])

def getBEDlen(BED):
	lengths = [x[2]-x[1] for x in BED]
	return sum(lengths)

# ---------------------------------------------------------------------------

# =============
# Get arguments
# =============

# Initiate the parser
parser = argparse.ArgumentParser()
parser.add_argument("-r", "--ref", help="reference genome assembly (multi fasta)", required=True)
parser.add_argument("-d", "--draft", help="draft genome assemblies (multi fasta)", required=True)
parser.add_argument("-o", "--output", help="Name of the outputed VCF", required=True)
parser.add_argument("-k", "--alignmentSizeKb", help="Minimum alignment size to have name correspondance (in kb)", type=int, default=20)
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
prefix=refName+"_"+draftName

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
	start=int(row[0])
	end=int(row[1])
	if start < end:
		BED2add = [[chr, start, end]]
	else:
		BED2add = [[chr, end, start]]		
	analysisMatrix[draftContigIndex][refChrIndex] = addBED(analysisMatrix[draftContigIndex][refChrIndex], BED2add)
coordsFile.close()
#os.remove(prefix+".coords")

# Convert all BED file in alignment length
for i in range(len(draftChr)):
	for j in range(len(refChr)):
		analysisMatrix[i][j] = getBEDlen(analysisMatrix[i][j])

# Get new contigs names
draftNewChr = []
print("--- Renamed contigs ---")
for i in range(len(draftChr)):
	correpondingNames=[]
	alignmentSizes=[]
	for j in range(len(refChr)):
		if analysisMatrix[i][j] >= alignmentSizeKb*1000:
			correpondingNames += [refChr[j]]
			alignmentSizes += [analysisMatrix[i][j]]

	if len(correpondingNames) == 0: # If no corresponding chromosome found
		newName=">"+draftChr[i]
		if args.alignmentSizeInfo:
			newName+=" length="+str(draftLen[i])
	else: 
		# Order chr according to decreasing alignment size
		correpondingNames=[correpondingNames[j] for j in decreasing_rank_simple(alignmentSizes)]
		alignmentSizes=[alignmentSizes[j] for j in decreasing_rank_simple(alignmentSizes)]

		newName=">"+draftChr[i]+"_"+"_".join(correpondingNames)
		if args.alignmentSizeInfo:
			newName+=" length="+str(draftLen[i])
			for j in range(len(correpondingNames)):
				newName+=" "+correpondingNames[j]+"="+str(alignmentSizes[j])
	print(draftChr[i]+"\t"+newName)
	draftNewChr += [newName]

# create new Fasta file with renamed contigs
out = open(outputPath, 'w')
for i in range(len(draftChr)):
	out.write(draftNewChr[i]+"\n")
	out.write(draftSeq[i]+"\n")
out.close()

print("New fasta written to: "+outputPath)
print("")

