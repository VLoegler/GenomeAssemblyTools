#!/usr/bin/env python
# -*- coding: utf8 -*-
#----------------------------------------------------------------------------
# Created By  : vloegler
# Created Date: 2022/09/15
# version ='1.0'
# ---------------------------------------------------------------------------
'''
This script check if all chromosomes of a reference are covered by
a draft assembly. 
It uses the nucmer (--maxmatch) and show-coords MUMmer4's functions

It takes as input :
	-r --ref: the reference genome assembly (multi fasta)
	-d --draft: draft genome assemblies (multi fasta, several allowed)
	-o --output: name of output tsv file
	-p --percentCovered: Percentage of each chromosome that have to be covered
	-m --mummerPath: Path to mummer function, if not in path
	-t --threads: Number of threads for nucmer

Output will be:
Assembly	NbChrCovered
Draft1		16
Draft2		15
Draft3		16
'''
# ---------------------------------------------------------------------------
import os
import sys
import argparse
import time
from datetime import datetime
from random import randint
from Tools import *
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
parser.add_argument("-d", "--draft", help="draft genome assemblies (multi fasta)", nargs='+', required=True)
parser.add_argument("-o", "--output", help="Output file", default = "")
parser.add_argument("-p", "--percentCovered", help="Percentage of each chromosome that have to be covered", type=int, default=80)
parser.add_argument("-m", "--mummerPath", help="Path to mummer function, if not in path", type=str, default="")
parser.add_argument("-t", "--threads", help="Number of threads for nucmer", type=int, default=1)


# Read arguments from the command line
args = parser.parse_args()

refPath=args.ref
draftPaths=args.draft
nbDraft=len(draftPaths)
outputPath=args.output
percent=int(args.percentCovered)
mummer=args.mummerPath
if mummer != "" and not mummer.endswith("/") :
	mummer += "/"
threads=args.threads

# Write header in output file
if outputPath != "":
	out = open(outputPath, 'w')
	out.write("Assembly\nNbChrPresent\n")
else:
	print("Assembly\tNbChrPresent")

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
# of chromosomes covered by the draft assembly
# ==================================================
for d in range(nbDraft):

	draftName=draftPaths[d].split("/")[-1]
	#print("\n\t--- Running for draft assembly "+draftName+" ---\n")
	start = time.time()

	# Align draft to ref
	prefix = "Alignment_" + refName + "_" + draftName + "_" + datetime.now().strftime("%Y_%m_%d_%H_%M_%S_%m_%f") + "_" + str(randint(0, 10000))
	getShowCoords(refPath, draftPaths[d], prefix, mummer, threads)
	end = time.time()
	#print("1/2 Alignment done: ran in "+str(round(end-start))+"s")

	start = time.time()
	# Get BED list of reference genome
	refBED = getBED(refPath)

	# ======================================
	# Get BED of alignments for each chromosome of the ref
	alignmentBEDs = [[] for i in range(len(refChr))]
	# Alignments BED will contain per Chromosome the BED of alignments of this chromosome on the draft genome
	coords = open(prefix+".coords", 'r')
	for line in coords.readlines():
		refStart = int(line.split("\t")[0])
		refEnd = int(line.split("\t")[1]) + 1
		Chr = line.split("\t")[7]
		ChrIndex = refChr.index(Chr)

		alignmentBEDs[ChrIndex] += [BEDcoordinates(id = Chr, start = refStart, end = refEnd)]
	coords.close()
	os.remove(prefix+".coords")

	# Convert all to BED
	alignmentBEDs = [BED(x) for x in alignmentBEDs]

	# Get coverage per chromosome
	refCoverage = []
	for i in range(len(refChr)):
		refCoverage += [refBED.getID(refChr[i]).overlapLen(alignmentBEDs[i], percent = True)]

	# Get number of chromosome covered at X%
	nbChrPresent = 0
	for i in range(len(refChr)):
		if refCoverage[i] >= percent:
			nbChrPresent += 1

	if outputPath != "":
		out.write(draftName + "\t" + str(nbChrPresent) + "\n")
	else: 
		print(draftName + "\t" + str(nbChrPresent))

if outputPath != "":
	out.close()




