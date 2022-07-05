#!/usr/bin/env python
# -*- coding: utf8 -*-
#----------------------------------------------------------------------------
# Created By  : vloegler
# Created Date: 2022/01/04
# version ='1.0'
# ---------------------------------------------------------------------------
'''
This script remove transloc haplotigs. In case the 2 haplotypes do not 
contain the same translocation, additional contig containing the translocation
can be in the draft assembly. Since this contig corresponds to a translocation,
it is fully covered by several other contigs of th assembly together. 

It takes as input :
	-d --draft: a draft genome assembly to reorder (multi fasta)
	-o --output: output prefix
	-b --blastPath: Path to blast+, if not in path


blast+ is used in this script. If blast+ is not in the path, the path can
be added to the variable blast. 
'''
# ---------------------------------------------------------------------------
import csv
import os
import sys
import argparse
from datetime import datetime
from random import randint
import re
import time
from Tools import *
# ---------------------------------------------------------------------------
# Definitions
def blastn(fasta1, fasta2, blastPath, out):
	blastCommand = blastPath + "blastn -query " + fasta1 + " -subject " + fasta2 + " -outfmt 6 -out " + out
	print(blastCommand)
	os.system(blastCommand)
# ---------------------------------------------------------------------------

# =============
# Get arguments
# =============

# Initiate the parser
parser = argparse.ArgumentParser()
parser.add_argument("-d", "--draft", help="draft genome assembly (multi fasta)", required=True)
parser.add_argument("-o", "--output", help="Prefix of the output file", required=True)
parser.add_argument("-b", "--blastPath", help="Path to blast+ function, if not in path", type=str, default="")

# Read arguments from the command line
args = parser.parse_args()

draftPath=args.draft
outputPath=args.output
blast=args.blastPath
if blast != "" and not blast.endswith("/") :
	blast += "/"

print("\n\t--- REMOVING TRANSLOC HAPLOTIGS ---\n")
print("Arguments detected:")
print("\t--draft:\t"+draftPath)
print("\t--output:\t"+outputPath)
print('')

# ===============
# Get Input files
# ===============
# Get draft contigs names and sequences
draftChr=[]
draftLine=[]
draftSeq=[]
seq=""
draft=open(draftPath, 'r')
for line in draft.readlines():
	if line.startswith(">"):
		draftLine += [line]
		draftChr += [line.split(">")[1].split(" ")[0].split("\t")[0].split("\n")[0]]
		if seq != "":
			draftSeq += [seq]
		seq=""
	else:
		seq += line.split("\n")[0]
draftSeq += [seq]
draft.close()
draftLen  = [len(x) for x in draftSeq]

blastResultsPath = datetime.now().strftime("%Y_%m_%d_%H_%M_%S_%m_%f") + "_" + str(randint(0, 10000)) + ".blastn"


# Run Blastn of draft against itself
start = time.time()
blastn(draftPath, draftPath, blast, blastResultsPath)
end = time.time()
print("Alignment done: ran in "+str(round(end-start))+"s")

# Get BED of the draft assemblies (BEDcoordinates without N nucleotides)
start = time.time()
draftBED = getBED(draftPath)
end = time.time()
print("Draft BED obtained: ran in "+str(round(end-start))+"s")

# Blast results analysis
# Blast results will be stored in the 2 level matrix analysisMatrix
# Level 1: Contig on which the contigs are aligned
# Level 2: BED of alignments covering contig1 by contig2
start = time.time()
x = len(draftChr)
y = len(draftChr)
analysisMatrix = []
for i in range(x):
	analysisMatrix.append([])
	for j in range(y):
		analysisMatrix[i].append([])

# Read blast out file
blastResultsFile = open(blastResultsPath)
blastResults = csv.reader(blastResultsFile, delimiter="\t")
for row in blastResults:
	contig1 = row[0]
	contig2 = row[1]
	if contig1 != contig2:
		index1 = draftChr.index(contig1)
		index2 = draftChr.index(contig2)
		startPos = int(row[6])
		endPos = int(row[7])+1
		coord2add = BEDcoordinates(id = contig1, start = startPos, end = endPos)
		analysisMatrix[index1][index2] += [coord2add]
blastResultsFile.close()
os.remove(blastResultsPath)

# Convert every list of BEDcoordinates to BED objects
for i in range(len(draftChr)):
	for j in range(len(draftChr)):
		analysisMatrix[i][j] = BED(analysisMatrix[i][j])
end = time.time()
print("Blast read: ran in "+str(round(end-start))+"s")

# Round 1: Compute coverage of each contig by all other contigs
start = time.time()
alignments = []
overlap = []
for i in range(len(draftChr)):
	# For each contig, the alignment list will contain the alignment BED of all other contigs on the contig
	contigs2align = list(range(len(draftChr)))
	contigs2align.remove(i)
	BED2sum = [analysisMatrix[i][x] for x in contigs2align]
	alignments += [BED(BED2sum)]

	# Convert each BED to overlap
	draftContigBED = draftBED.getID(draftChr[i]) # Get BED of contig
	overlap += [draftContigBED.overlapLen(alignments[i], percent = True)] # Get overlap percentage of alignment on contig

# While there is contigs covered on more than 90%, remove contig and recompute the whole process
contigsRemoved = []
contigsRemovedCoverage = []
while max(overlap) >= 95:
	contigsRemoved += [overlap.index(max(overlap))]
	contigsRemovedCoverage += [max(overlap)]

	alignments = [[] for i in range(len(draftChr))]
	overlap = [0 for i in range(len(draftChr))]
	# For each contig, the alignment list will contain the alignment BED of all other contigs on the contig
	for i in range(len(draftChr)):
		if i not in contigsRemoved:
			# For each contig, the alignment list will contain the alignment BED of all other contigs on the contig
			contigs2align = list(range(len(draftChr)))
			contigs2align.remove(i)
			for j in contigsRemoved:
				contigs2align.remove(j)
			BED2sum = [analysisMatrix[i][x] for x in contigs2align]
			alignments[i] = BED(BED2sum)

			# Convert each BED to overlap
			draftContigBED = draftBED.getID(draftChr[i])
			overlap[i] = draftContigBED.overlapLen(alignments[i], percent = True)

end = time.time()
print("Coverage analysis completed: ran in "+str(round(end-start))+"s")

# Write output to 2 files: PREFIX.NoTranslocHaplotigs.fasta and PREFIX.removedTranslocHaplotigs.fasta
output1=open(outputPath+".NoTranslocHaplotigs.fasta", 'w')
if len(contigsRemoved) > 0 :
	output2 = open(outputPath+".removedTranslocHaplotigs.fasta", "w")
for i in range(len(draftChr)):
	seq=re.sub("(.{80})", "\\1\n", draftSeq[i], 0, re.DOTALL)
	if i in contigsRemoved:
		output2.write(draftLine[i])
		output2.write(seq+"\n")
	else:
		output1.write(draftLine[i])
		output1.write(seq+"\n")
output1.close()
if len(contigsRemoved) > 0 :
	output2.close()
	print("\nContigs removed: ")
	print("\tContig\t% covered")
	for i in range(len(contigsRemoved)):
		print("\t" + draftChr[contigsRemoved[i]] + "\t" + str(contigsRemovedCoverage[i]))
else:
	print("\nNo contig removed")
