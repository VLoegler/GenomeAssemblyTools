#!/usr/bin/env python
# -*- coding: utf8 -*-
#----------------------------------------------------------------------------
# Created By  : vloegler
# Created Date: 2022/03/04
# version ='1.0'
# ---------------------------------------------------------------------------
'''
This script check looks for wrong translocation in several assemblies made
from the the same sequencing data. A translocation is considered as one
contig aligning on more than X kb on two ore more reference chromosomes. 
If another assembly have a contig aligning on 50% on both Transloc and non-
traslocated part of the reference chromosome, the translocation is considered
as wrong. 
It's better to take a reference with telomere, TEs and repeats masked. 

Translocation detected :
		_________________________
		|          |           /|
		|          |          / |
		|          |         /  |
		|          |        /   |
contig	|     /    |            |
		|    /     |            |
		|   /      |            |
		|  /       |            |
		| /        |            |
		|/_________|____________|
		  Ref CHR1	Ref CHR2

Transloc wrong if  :
			____________
			|         /|
			|        / |
			|       /  |
			|      /   |
contig		|     /    |
of			|    /     |
another		|   /      |
assembly	|  /       |
			| /        |
			|/_________|
		  Ref CHR1 or CHR2


It takes as input :
	-r --ref: the reference genome assembly (multi fasta)
	-d --draft: draft genome assemblies from the same strain of sequencing dataset (multi fasta, several allowed)
	-o --output: name of output tsv file
	-k --alignmentSizeKb: Minimum alignment size to have name correspondance (in kb)
	-m --mummerPath: Path to mummer function, if not in path
	-t --threads: Number of threads for nucmer

Output will be:
Assembly	WrongTransloc
Draft1		False
Draft2		False
Draft3		True
'''
# ---------------------------------------------------------------------------
from Tools import *
import os
import sys
import argparse
import time
import csv
import copy
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

# ---------------------------------------------------------------------------

# =============
# Get arguments
# =============

# Initiate the parser
parser = argparse.ArgumentParser()
parser.add_argument("-r", "--ref", help="reference genome assembly (multi fasta)", required=True)
parser.add_argument("-d", "--draft", help="draft genome assemblies from the same strain of sequencing dataset (multi fasta)", nargs='+', required=True)
parser.add_argument("-o", "--output", help="Prefix of the output", required=True)
parser.add_argument("-k", "--alignmentSizeKb", help="Minimum alignment size to detect translocation (in kb)", type=int, default=30)
parser.add_argument("-m", "--mummerPath", help="Path to mummer function, if not in path", type=str, default="")
parser.add_argument("-t", "--threads", help="Number of threads for nucmer", type=int, default=1)


# Read arguments from the command line
args = parser.parse_args()

refPath=args.ref
draftPaths=args.draft
nbDraft=len(draftPaths)
outputPath=args.output
alignmentSizeKb=args.alignmentSizeKb
mummer=args.mummerPath
if mummer != "" and not mummer.endswith("/") :
	mummer += "/"
threads=args.threads

# Write header in output file
out = open(outputPath+".tsv", 'w')
out.write("Assembly\tWrongTransloc\n")
outDetailed = open(outputPath+".detailed.tsv", 'w')
outDetailed.write("Assembly\tWrongTransloc\tContig\tContigSize\tNbRefChrAligned\tRefChr1\tAlignmentSize1\tRefChr2\tAlignmentSize2\n")

# ========================================
# Get reference chromosome name and length
# ========================================
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

# Get BED list of reference genome
print("\nGetting non N positions of reference in BED")
startTime = time.time()
refBED = getBED(refPath)
endTime = time.time()
print("\tObtained in "+str(round(endTime-startTime))+"s")

# ==================================
# Align each draft against reference
# ==================================
print("\nRunning alignments for draft assemblies")
startTime = time.time()
for d in range(nbDraft):
	draftName=draftPaths[d].split("/")[-1]
	# Align draft to ref
	prefix = draftName+"_VS_"+refName
	getShowCoords(refPath, draftPaths[d], prefix, mummer, threads)
endTime = time.time()
print("\tAlignments ran in "+str(round(endTime-startTime))+"s")


# ============================================
# For each assembly, search for translocations
# ============================================
startTime = time.time()
for d in range(nbDraft):

	draftName=draftPaths[d].split("/")[-1]
	print("Searching for translocations in " + draftName)

	# Assuming at first that there is no wrong translocation
	wrongTransloc = False

	# Get draft contig name
	draftChr=[]
	draftSeq=[]
	seq=""
	draft=open(draftPaths[d], 'r')
	for line in draft.readlines():
		if line.startswith(">"):
			draftChr += [line.strip().split(">")[1].split(" ")[0].split("\t")[0]]
			if seq != "":
				draftSeq += [seq]
			seq = ""
		else:
			seq += line.strip()
	draftSeq += [seq]
	draft.close()
	draftLen=[len(x) for x in draftSeq]

	# Get prefix of alignment file
	prefix = draftName+"_VS_"+refName

	# Create a matrix which will contain BED objects of aligned regions
	# First dimension: draft contigs
	# Second dimension: reference chromosomes
	x = len(draftChr)
	y = len(refChr)
	analysisMatrix = []
	for i in range(x):
		analysisMatrix.append([])
		for j in range(y):
			analysisMatrix[i].append([])

	coordsFile = open(prefix+".coords")
	coords = csv.reader(coordsFile, delimiter="\t")
	for row in coords:
		draftContigIndex=draftChr.index(row[8])
		refChrIndex=refChr.index(row[7])
		# get BED to add
		chr=row[7]
		startPos=int(row[0])
		endPos=int(row[1])
		coord2add=BEDcoordinates(id = chr, start = startPos, end = endPos)
		# Substract BED
		analysisMatrix[draftContigIndex][refChrIndex] += [coord2add]
	coordsFile.close()

	# Convert all to BED objects
	for i in range(len(draftChr)):
		for j in range(len(refChr)):
			analysisMatrix[i][j] = BED(analysisMatrix[i][j])

	# For each point of the matrix, get the size of aligned regions
	lenMatrix = copy.deepcopy(analysisMatrix)
	for i in range(len(draftChr)):
		for j in range(len(refChr)):
			# Get length of each BED
			lenMatrix[i][j] = lenMatrix[i][j].getLen()

	# Run over draft contigs to see if a contig aligns on on several reference chromosomes
	for i in range(len(draftChr)):
		chrAboveX = [lenMatrix[i].index(x) for x in lenMatrix[i] if x >= alignmentSizeKb*1e3]
		# If several contigs are aligned, that is a transloc
		if len(chrAboveX) > 1:
			# Get the BED coordinates of aligned part of chromosome and non aligned part
			refChrIndices = chrAboveX
			aligned = [analysisMatrix[i][j] for j in chrAboveX]
			nonaligned = []

			for j in range(len(chrAboveX)):
				# Get reference BED coordinates for the chromosome
				chrIndex = chrAboveX[j]
				chrBED = refBED.getID(refChr[chrIndex])
				# Substract aligned part
				chrBED.substractBED(aligned[j])
				nonaligned += [chrBED.copy()]

			# =================================================================================
			# Check if duplication: 
			# If another contig of the same assembly reject the translocation (the contig
			# aligns both on aligned and non aligned part), the event is probably a duplication
			# of a region. Since the translocation is spuriously rejected, the transloc is
			# not considered as "wrong". 
			# =================================================================================
			duplication = False

			# For each contig, check if the contig aligns on aligned part of the translocation
			# In that case, it is a duplication event
			for i2 in range(len(draftChr)):
				if i != i2 :
					for j in range(len(chrAboveX)):
						refChrIndex = refChrIndices[j]
						contigBED = analysisMatrix[i2][refChrIndex]

						# Si BED overlap + de 95% de aligned : DUPLICATION !
						if aligned[j].overlapLen(contigBED, percent = True) >= 95:
							duplication = True
							print("Translocation detected due to a duplication event")
							print("\t" + draftChr[i] + " aligns on")
							for j in chrAboveX:
								print("\t\t" + refChr[j] + " (alignment on " + str(round(lenMatrix[i][j] / 1000)) + "kb)")

			if duplication == False:
				# ===================================================================
				# check if a contig of another assembly aligns on 95% of both aligned 
				# and non aligned part of the ref chromosome. If it is the case, the 
				# translocation is considered as wrong
				# ===================================================================
				# For file in all assemblies: 
				for d2 in range(nbDraft):
					if d != d2: # If assembly to compare is not the assembly with transloc found
						draftName2=draftPaths[d2].split("/")[-1]
						prefix2 = draftName2+"_VS_"+refName

						# Get draft contig name
						draftChr2=[]
						draft2=open(draftPaths[d2], 'r')
						for line in draft2.readlines():
							if line.startswith(">"):
								draftChr2 += [line.strip().split(">")[1].split(" ")[0].split("\t")[0]]
						draft2.close()

						# Create a matrix which will contain BED objects of aligned regions
						# First dimension: draft contigs
						# Second dimension: reference chromosomes
						x2 = len(draftChr2)
						y2 = len(refChr)
						analysisMatrix2 = []
						for i2 in range(x2):
							analysisMatrix2.append([])
							for j2 in range(y2):
								analysisMatrix2[i2].append([])
								
						coordsFile = open(prefix2+".coords")
						coords = csv.reader(coordsFile, delimiter="\t")
						for row in coords:
							draftContigIndex2=draftChr2.index(row[8])
							refChrIndex=refChr.index(row[7])
							# get BED to add
							chr=row[7]
							start=int(row[0])
							end=int(row[1])
							coord2add=BEDcoordinates(id = chr, start = start, end = end)
							# Substract BED
							analysisMatrix2[draftContigIndex2][refChrIndex] += [coord2add]
						coordsFile.close()

						# Convert all to BED objects
						for i2 in range(len(draftChr2)):
							for j2 in range(len(refChr)):
								analysisMatrix2[i2][j2] = BED(analysisMatrix2[i2][j2])

						# For each contig, check if the contig aligns on both aligned and non align part of the translocation
						for i2 in range(len(draftChr2)):
							for j in range(len(chrAboveX)):
								refChrIndex = refChrIndices[j]
								contigBED = analysisMatrix2[i2][refChrIndex]
								#print(draftChr2[i2] + "\t" + refChr[refChrIndex] + "\tAligned: " + str(aligned[j].overlapLen(contigBED, percent = True)) + "\tNonAligned: " + str(nonaligned[j].overlapLen(contigBED, percent = True)))

								# Si BED overlap + de 95% de aligned et non aligned : WRONG TRANSLOC !
								if aligned[j].overlapLen(contigBED, percent = True) >= 95 and nonaligned[j].overlapLen(contigBED, percent = True) >= 95:
									if not wrongTransloc:
										print("\nWrong translocation found:")
										print("\t" + draftChr[i] + " aligns on")
										for k in chrAboveX:
											print("\t\t" + refChr[k] + " (alignment on " + str(round(lenMatrix[i][k] / 1000)) + "kb)")
										print("Proved wrong by draft" + draftName2)
										print("draftChr\trefChr\tAlignedSiteCoverage\tNonAlignedSiteCoverage")
										print(draftChr2[i2] + "\t" + refChr[refChrIndex] + "\t" + str(round(aligned[j].overlapLen(contigBED, percent = True), 2)) + "\t" + str(round(nonaligned[j].overlapLen(contigBED, percent = True), 2)) + "\n")

										# Write transloc in detailed output file
										index1 = chrAboveX[0]
										index2 = chrAboveX[1]
										outDetailed.write(draftName + "\t" + str(wrongTransloc) + "\t" + draftChr[i] + "\t" + str(draftLen[i]) + "\t" + str(len(chrAboveX)) + "\t" + refChr[index1] + "\t" + str(lenMatrix[i][index1]) + "\t" + refChr[index2] + "\t" + str(lenMatrix[i][index2]) + "\n")
									wrongTransloc = True



	# Write presence or not of a wrong translocation to the output file
	out.write(draftName + "\t" + str(wrongTransloc) + "\n")
	# If no transloc found, write in detailed output file
	if not wrongTransloc:
		outDetailed.write(draftName + "\t" + str(wrongTransloc) + "\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n")
out.close()
outDetailed.close()
endTime = time.time()
print("\nAnalysis ran in "+str(round(endTime-startTime))+"s")


# ==================
# Remove coord files
# ==================
for d in range(nbDraft):
	draftName=draftPaths[d].split("/")[-1]
	prefix = draftName+"_VS_"+refName
	os.remove(prefix+".coords")

