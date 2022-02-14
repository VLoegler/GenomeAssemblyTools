#!/usr/bin/env python
# -*- coding: utf8 -*-
#----------------------------------------------------------------------------
# Created By  : vloegler
# Created Date: 2022/02/11
# version ='1.0'
# ---------------------------------------------------------------------------
'''
This script check how many contigs are needed to cover X% of the reference

It takes as input :
	-r --ref: the reference genome assembly (multi fasta)
	-d --draft: draft genome assemblies (multi fasta, several allowed)
	-o --output: name of output tsv file
	-p --percentCovered: Percentage of the genome which has to be covered
	-m --mummerPath: Path to mummer function, if not in path
	-e --enableMultiprocessing: Enable multiprocessing, highly speed up the script

Output will be:
Assembly	NbContigsToX
Draft1		1
Draft2		1
Draft3		2
'''
# ---------------------------------------------------------------------------
import csv
import os
import sys
import argparse
from datetime import datetime
from random import randint
# ---------------------------------------------------------------------------
# Definitions

def rank_simple(vector):
	return sorted(range(len(vector)), key=vector.__getitem__)

def decreasing_rank_simple(vector):
	return sorted(range(len(vector)), key=vector.__getitem__)[::-1]

def createFasta(name: str, IDs: list, sequences: list):
	if len(IDs) != len(sequences) :
		raise Error("List of IDs and sequences must have the same length. Fasta file connt be created.")
	fasta = open(name, "w")
	for i in range(len(IDs)):
		fasta.write(">"+IDs[i]+"\n")
		fasta.write(sequences[i]+"\n")
	fasta.close()

def getPercentageCovered(refPath, draftPath, mummerPath):
	# Run dnadiff command
	print("Running dnadiff between:")
	print(refPath)
	print(draftPath)
	prefix = datetime.now().strftime("%Y_%m_%d_%H_%M_%S_%m_%f") + "_dnadiff_" + str(randint(0, 10000))
	dnadiffShellCommand=mummerPath+"dnadiff --prefix "+prefix+" "+refPath+" "+draftPath
	os.system(dnadiffShellCommand)
	# Read report file to get the percentage of base aligned
	report = open(prefix + ".report", "r")
	for line in report.readlines():
		if line.startswith("AlignedBases"):
			percentCovered = float(line.split("%")[0].split("(")[1])
			break
	# Remove dnadiff files
	os.system("rm -f "+prefix+".*")

	return percentCovered



# ---------------------------------------------------------------------------

# =============
# Get arguments
# =============

# Initiate the parser
parser = argparse.ArgumentParser()
parser.add_argument("-r", "--ref", help="reference genome assembly (multi fasta)", required=True)
parser.add_argument("-d", "--draft", help="draft genome assemblies (multi fasta)", nargs='+', required=True)
parser.add_argument("-o", "--output", help="Name of the outputed VCF", required=True)
parser.add_argument("-p", "--percentCovered", help="Percentage of the genome which has to be covered", type=int, default=99)
parser.add_argument("-m", "--mummerPath", help="Path to mummer function, if not in path", type=str, default="")
parser.add_argument("-e", "--enableMultiprocessing", help="Enable multiprocessing, highly speed up the script", type=bool, default=True)


# Read arguments from the command line
args = parser.parse_args()

refPath=args.ref
draftPaths=args.draft
nbDraft=len(draftPaths)
outputPath=args.output
percent=args.percentCovered
mummer=args.mummerPath
if mummer != "" and not mummer.endswith("/") :
	mummer += "/"
enableMultiprocessing=args.enableMultiprocessing
if enableMultiprocessing:
	import multiprocessing as mp

# Write header in output file
out = open(outputPath, 'w')
out.write("Assembly\tNbContigsTo" + str(percent) + "\n")

# Iterate over each draft assembly
for d in range(nbDraft):
	draftName = draftPaths[d].split("/")[-1]

	# =====================================
	# Get draft contigs names and sequences
	# =====================================
	draftChr=[]
	draftLine=[]
	draftSeq=[]
	seq=""
	draft=open(draftPaths[d], 'r')
	for line in draft.readlines():
		if line.startswith(">"):
			draftLine += [line]
			draftChr += [line.strip().split(">")[1].split(" ")[0].split("\t")[0]]
			if seq != "":
				draftSeq += [seq]
			seq=""
		else:
			seq += line.strip()
	draftSeq += [seq]
	draft.close()


	# ===============================================
	# Using MUMmer4 to align each contig to reference
	# ===============================================

	# Get the fraction of the reference covered for each contig
	toRun = [] # This list will contain the arguments required to run the getPercentageCovered function
	for i in range(len(draftChr)):

		# Create a fasta file with a single contig
		fastaName = draftName+"_"+draftChr[i]+".fasta"
		createFasta(fastaName, [draftChr[i]], [draftSeq[i]])

		# Use MUMmer4's dnadiff to get the percentage of the reference covered by the contig
		toRun += [(refPath, fastaName, mummer)] 

	# Run dnadiff to get cevorage of the reference
	if enableMultiprocessing:
		pool = mp.Pool(mp.cpu_count())
		draftCoverage = pool.starmap(getPercentageCovered, toRun)
	else:
		draftCoverage = []
		for run in toRun:
			draftCoverage += [getPercentageCovered(run[0], run[1], run[2])]
	print(draftCoverage)

	# Remove fasta files
	for run in toRun:
		fasta = run[1]
		os.remove(fasta)

	# Order contigs base on decreasing prcentage of reference covered
	order = decreasing_rank_simple(draftCoverage)

	# ================================================================
	# Get the number of contigs necessary to cover X% of the reference
	# ================================================================

	# First contigs
	n = 1
	index = order[0]
	ID = [draftChr[index]]
	sequence = [draftSeq[index]]
	fastaName = draftName+"_"+str(n)+"_contigs.fasta"
	createFasta(fastaName, ID, sequence)
	toRun = [(refPath, fastaName, mummer)]

	# Sequentially add all contigs
	for i in range(1,len(draftChr)):
		n += 1
		index = order[0:n]
		# Create Fasta File with contigs 1 to n
		IDs = []
		sequences = []
		for i in index:
			IDs += [draftChr[i]]
			sequences += [draftSeq[i]]
		fastaName = draftName+"_"+str(n)+"_contigs.fasta"
		createFasta(fastaName, IDs, sequences)
		toRun += [(refPath, fastaName, mummer)] 

	# Run dnadiff
	if enableMultiprocessing:
		pool = mp.Pool(mp.cpu_count())
		increasingDraftCoverage = pool.starmap(getPercentageCovered, toRun)
	else:
		increasingDraftCoverage = []
		for run in toRun:
			increasingDraftCoverage += [getPercentageCovered(run[0], run[1], run[2])]
	print(increasingDraftCoverage)

	# Remove fasta files
	for run in toRun:
		fasta = run[1]
		os.remove(fasta)

	# Get the number of contigs to reach X%
	referenceCovered = 0
	nbContigs = 0
	for c in increasingDraftCoverage:
		nbContigs += 1
		if c >= percent:
			break
	if c >= percent:
		n = nbContigs
	else:
		n = "NotReached"

	out.write(draftName + "\t" + str(n) + "\n")

out.close()


