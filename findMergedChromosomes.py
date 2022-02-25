#!/usr/bin/env python
# -*- coding: utf8 -*-
#----------------------------------------------------------------------------
# Created By  : vloegler
# Created Date: 2022/02/24
# version ='1.0'
# ---------------------------------------------------------------------------
'''
This script check if there are contigs corresponding to merged chromosomes
of a reference. It looks if a single contigs aligns on more than 90% of 2 
reference chromosomes. 
It uses the nucmer (--maxmatch) and show-coords MUMmer4's functions

It takes as input :
	-r --ref: the reference genome assembly (multi fasta)
	-d --draft: draft genome assemblies (multi fasta, several allowed)
	-o --output: name of output tsv file
	-p --percentToMatch: Percentage a contig has to cover on each chromosome to be considered as a merge
	-m --mummerPath: Path to mummer function, if not in path
	-t --threads: Number of threads for nucmer

Output will be:
Assembly	NbMergedContigs
Draft1		0
Draft2		0
Draft3		1
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
	nucmerShellCommand=mummerPath+"nucmer -t "+str(threads)+" --maxmatch --prefix "+prefix+" "+refPath+" "+draftPath
	os.system(nucmerShellCommand)
	# Run show coords
	showcoordsShellCommand=mummerPath+"show-coords -TH "+prefix+".delta > "+prefix+".coords"
	os.system(showcoordsShellCommand)
	# Remove nucmer delta file
	os.remove(prefix+".delta")

def getBed(fastaPath):
	# Get reference chromosome name and length
	Chr=[]
	Seq=[]
	seq=""
	fasta=open(fastaPath, 'r')
	for line in fasta.readlines():
		if line.startswith(">"):
			Chr += [line.strip().split(">")[1].split(" ")[0].split("\t")[0]]
			if seq != "":
				Seq += [seq]
			seq=""
		else:
			seq += line.strip()
	Seq += [seq]
	fasta.close()
	Len=[len(x) for x in Seq]
	# Create the BED file for each chromosome [ChrID, Start, End]
	# Start is included, End is excluded
	BED = []
	for i in range(len(Chr)):
		BED += [[Chr[i], 1, Len[i]+1]]
	return BED

def substractBEDcoordinates(b1, b2):
	# Get info from b1
	chr1=b1[0]
	start1=b1[1]
	end1=b1[2]
	# Get info from b2
	chr2=b2[0]
	start2=b2[1]
	end2=b2[2]
	# Store b1-b2 in newb BED object
	newb=[]
	if chr2 == chr1:
		if start2 <= start1:
			if end2 <= start1:
				newb = [b1]
			elif end2 < end1:
				newb = [[chr1, end2, end1]]
			else: 
				newb = [] # the BED coordinates are deleted
		elif start2 < end1:
			if end2 < end1:
				newb = [[chr1, start1, start2], [chr1, end2, end1]]
			else: 
				newb = [[refChr, start1, start2]]
		else:
			newb = [b1]
	else: 
		newb = [b1]
	# Return the substraction as BED object (list of BED coordinates)
	return newb

def substractBED(BED1, BED2):
	# Remove each coordinates of BED2 to all BED1 coordinates
	for b2 in BED2:
		chr2=b2[0]
		start2=b2[1]
		end2=b2[2]
		# newBED will store the temporary subtracted bed
		newBED = []
		# Substract b2 to all b1 coodinates
		for b1 in BED1:
			chr1=b1[0]
			start1=b1[1]
			end1=b1[2]
			newb1=substractBEDcoordinates(b1, b2)
			newBED += newb1
		# Replace old BED1 by the new one
		BED1 = newBED.copy()
	return BED1

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
parser.add_argument("-d", "--draft", help="draft genome assemblies (multi fasta)", nargs='+', required=True)
parser.add_argument("-o", "--output", help="Name of the outputed VCF", required=True)
parser.add_argument("-p", "--percentToMatch", help="Percentage a contig has to cover on each chromosome to be considered as a merge", type=int, default=90)
parser.add_argument("-m", "--mummerPath", help="Path to mummer function, if not in path", type=str, default="")
parser.add_argument("-t", "--threads", help="Number of threads for nucmer", type=int, default=1)


# Read arguments from the command line
args = parser.parse_args()

refPath=args.ref
draftPaths=args.draft
nbDraft=len(draftPaths)
outputPath=args.output
percent=args.percentToMatch
mummer=args.mummerPath
if mummer != "" and not mummer.endswith("/") :
	mummer += "/"
threads=args.threads

# Write header in output file
out = open(outputPath, 'w')
out.write("Assembly\tNbMergedContigs\n")

# ========================================
# Get reference chromosome name and length
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

refName=refPath.split("/")[-1]

# ==================================================
# Iterate over each draft assembly to get the number 
# of contigs that cover X% of the reference
# ==================================================
for d in range(nbDraft):

	draftName=draftPaths[d].split("/")[-1]
	print("\n\t--- Running for draft assembly "+draftName+" ---\n")
	startTime = time.time()

	# Get draft contig name
	draftChr=[]
	draft=open(draftPaths[d], 'r')
	for line in draft.readlines():
		if line.startswith(">"):
			draftChr += [line.strip().split(">")[1].split(" ")[0].split("\t")[0]]
	draft.close()

	# Align draft to ref
	prefix = draftName+"_VS_"+refName
	getShowCoords(refPath, draftPaths[d], prefix, mummer, threads)
	endTime = time.time()
	print("1/2 Alignment done: ran in "+str(round(endTime-startTime))+"s")

	startTime = time.time()
	# Get BED list of reference genome
	refBED = getBed(refPath)

	# Create a matrix which will contain BED objects of uncovered regions
	# First dimension: draft contigs
	# Second dimension: reference chromosomes
	x = len(draftChr)
	y = len(refChr)
	analysisMatrix = []
	for i in range(x):
		analysisMatrix.append([])
		for j in range(y):
			analysisMatrix[i].append([refBED[j]])

	# For each line in the alignment coords file, substract BED files
	coordsFile = open(prefix+".coords")
	coords = csv.reader(coordsFile, delimiter="\t")
	for row in coords:
		draftContigIndex=draftChr.index(row[8])
		refChrIndex=refChr.index(row[7])
		# get BED to substract
		chr=row[7]
		start=int(row[0])
		end=int(row[1])
		BED2substract=[[chr, start, end]]
		# Substract BED
		analysisMatrix[draftContigIndex][refChrIndex] = substractBED(analysisMatrix[draftContigIndex][refChrIndex], BED2substract)
	coordsFile.close()
	os.remove(prefix+".coords")

	# For each point of the matrix, get the percentage of reference chromosome covered
	for i in range(len(draftChr)):
		for j in range(len(refChr)):
			# Get uncovered length
			uncoveredLen = getBEDlen(analysisMatrix[i][j])
			# Get reference chromosome size
			refChrLen = refLen[j]
			# Compute proportion of chromosome covered
			analysisMatrix[i][j] = 100 * ( refChrLen - uncoveredLen ) / refChrLen

	# For each draft contig, check if it is not 2 merged reference chromosomes
	nbContigsMerged = 0
	for i in range(len(draftChr)):
		# If coverage is superior to 90% for more than one ref chromosome: it is a merged contig
		supToX = len([x for x in analysisMatrix[i] if x >= percent])
		if supToX > 1:
			nbContigsMerged += 1

	endTime = time.time()
	print("2/2 Coordinates analsis done: ran in "+str(round(endTime-startTime))+"s\n")

	out.write(draftName + "\t" + str(nbContigsMerged) + "\n")
	print("Assembly\tNbMergedContigs")
	print(draftName + "\t" + str(nbContigsMerged))

out.close()
