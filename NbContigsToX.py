#!/usr/bin/env python
# -*- coding: utf8 -*-
#----------------------------------------------------------------------------
# Created By  : vloegler
# Created Date: 2022/02/11
# version ='1.0'
# ---------------------------------------------------------------------------
'''
This script check how many contigs are needed to cover X% of the reference
It uses the nucmer (--maxmatch) and show-coords MUMmer4's functions

It takes as input :
	-r --ref: the reference genome assembly (multi fasta)
	-d --draft: draft genome assemblies (multi fasta, several allowed)
	-o --output: name of output tsv file
	-p --percentCovered: Percentage of the genome which has to be covered
	-m --mummerPath: Path to mummer function, if not in path
	-t --threads: Number of threads for nucmer

Output will be:
Assembly	NbContigsToX
Draft1		16
Draft2		17
Draft3		16
'''
# ---------------------------------------------------------------------------
import os
import sys
import argparse
import time
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

def getBEDlen(BED):
	lengths = [x[2]-x[1] for x in BED]
	return sum(lengths)

def getCoverage(refBED, showCoordsPath, contigList:list):
	BED = refBED.copy()
	coords = open(showCoordsPath, 'r')
	for line in coords.readlines():
		refStart = int(line.split("\t")[0])
		refEnd = int(line.split("\t")[1]) + 1
		refChr = line.split("\t")[7]
		draftChr = line.strip().split("\t")[8]

		if draftChr in contigList:
			newBED = []
			for b in BED:
				if b[0] == refChr:
					if refStart <= b[1]:
						if refEnd <= b[1]:
							newBED += [b]
						elif refEnd < b[2]:
							newBED += [[refChr, refEnd, b[2]]]
						else: 
							pass # the BED coordinates are deleted
					elif refStart < b[2]:
						if refEnd < b[2]:
							newBED += [[refChr, b[1], refStart]]
							newBED += [[refChr, refEnd, b[2]]]
						else: 
							newBED += [[refChr, b[1], refStart]]
					else:
						newBED += [b]
				else: newBED += [b]
			BED = newBED.copy()
	coords.close()

	# Get total length of ref Genome
	totalLen = sum(refLen)
	uncoveredLen = getBEDlen(BED)
	coveredLen = totalLen - uncoveredLen
	percent = 100 * coveredLen / totalLen
	return(percent)
# ---------------------------------------------------------------------------

# =============
# Get arguments
# =============

# Initiate the parser
parser = argparse.ArgumentParser()
parser.add_argument("-r", "--ref", help="reference genome assembly (multi fasta)", required=True)
parser.add_argument("-d", "--draft", help="draft genome assemblies (multi fasta)", nargs='+', required=True)
parser.add_argument("-o", "--output", help="Name of the outputed VCF", required=True)
parser.add_argument("-p", "--percentCovered", help="Percentage of the genome which has to be covered", type=int, default=95)
parser.add_argument("-m", "--mummerPath", help="Path to mummer function, if not in path", type=str, default="")
parser.add_argument("-t", "--threads", help="Number of threads for nucmer", type=int, default=1)


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
threads=args.threads

# Write header in output file
out = open(outputPath, 'w')
out.write("Assembly\tNbContigsTo" + str(percent) + "\n")

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
	start = time.time()

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
	end = time.time()
	print("1/2 Alignment done: ran in "+str(round(end-start))+"s")


	start = time.time()
	# Get BED list of reference genome
	refBED = getBed(refPath)

	# ======================================
	# Get reference coverage for each contig
	toRun = []
	for contig in draftChr:
		toRun += [[refBED, prefix+".coords", [contig]]]

	# Get coverage per contig
	draftCoverage = []
	for run in toRun:
		draftCoverage += [getCoverage(run[0], run[1], run[2])]

	# Order contigs by decreasing coverage of reference assembly
	order = decreasing_rank_simple(draftCoverage)

	# =========================================================
	# Get reference coverage by adding sequentially each contig
	toRun = []
	# First contig
	index = order[0]
	contig = [draftChr[index]]
	toRun += [[refBED, prefix+".coords", contig]]

	# Sequentially add all contigs
	for i in range(1,len(draftChr)):
		index = order[0:i+1]
		contigs = [draftChr[j] for j in index]
		toRun += [[refBED, prefix+".coords", contigs]]

	# Get coverage
	increasingDraftCoverage = []
	for run in toRun:
		increasingDraftCoverage += [getCoverage(run[0], run[1], run[2])]

	nbContigs = 0
	for c in increasingDraftCoverage:
		nbContigs += 1
		if c >= percent:
			break
	if c >= percent:
		n = nbContigs
	else:
		n = "NotReached"
	os.remove(prefix+".coords")

	end = time.time()
	print("2/2 Coordinates analsis done: ran in "+str(round(end-start))+"s\n")

	out.write(draftName + "\t" + str(n) + "\n")
	print("Assembly\tNbContigsTo" + str(percent))
	print(draftName + "\t" + str(n))

out.close()




