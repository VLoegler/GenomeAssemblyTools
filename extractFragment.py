#!/usr/bin/env python
# -*- coding: utf8 -*-
#----------------------------------------------------------------------------
# Created By  : vloegler
# Created Date: 2022/09/15
# version ='1.2'
# ---------------------------------------------------------------------------
'''
This script extract a sequence from a fasta file

It takes as input :
	-f --fasta: fasta file
	-id --seqid: ID of the sequence, or multiple IDs separated by space
	-s --start: Start position (inclusive, 1-based)
	-e --end: End position (inclusive, 1-based)
	-o --output: output file
'''
# ---------------------------------------------------------------------------
import os
import sys
import argparse
import re
# ---------------------------------------------------------------------------

# =============
# Get arguments
# =============

# Initiate the parser
parser = argparse.ArgumentParser()
parser.add_argument("-f", "--fasta", help="fasta file", required=True)
parser.add_argument("-o", "--output", help="output file", default = "")
parser.add_argument("-id", "--seqid", help="ID of the sequence", nargs='+', required=True)
parser.add_argument("-s", "--start", help="Start position (inclusive, 1-based)", type = int, default = 0)
parser.add_argument("-e", "--end", help="End position (inclusive, 1-based)", type = int, default = 0)

# Read arguments from the command line
args = parser.parse_args()

fastaPath = args.fasta
outputPath=args.output
seqidList=args.seqid
start=int(args.start)
end=int(args.end)

# ===============
# Get Input files
# ===============

# Get reference chromosomes names & sequences
Chr=[]
Seq=[]
seq=""
fasta=open(fastaPath, 'r')
for line in fasta.readlines():
	if line.startswith(">"):
		Chr += [line.split(">")[1].split(" ")[0].strip()]
		if seq != "":
			Seq += [seq]
		seq=""
	else:
		seq += line.split("\n")[0]
Seq += [seq]
fasta.close()

# Open output file if necessary
if outputPath != "":
	out = open(outputPath, "w")

for seqid in seqidList:
	if seqid not in Chr:
		raise ValueError("Sequence ID is not in fasta file")

	targetSeq = Seq[Chr.index(seqid)]

	if start != 0 or end != 0:
		if start < 1 or start > len(targetSeq) or end < 1 or end > len(targetSeq) or end <= start:
			raise ValueError("Bad coordinates. Start and End must be 1-based, not higher than sequence length and starts < end.")

	if start != 0 or end != 0:
		extractedSeq = targetSeq[(start-1):end]
	else:
		extractedSeq = targetSeq

	if outputPath != "":
		if start != 0 or end != 0:
			out.write(">" + seqid + "_" + str(start) + ".." + str(end) + "\n")
		else:
			out.write(">" + seqid + "\n")
		seq=re.sub("(.{80})", "\\1\n", extractedSeq, 0, re.DOTALL)
		out.write(seq+"\n")
	else:
		if start != 0 or end != 0:
			print(">" + seqid + "_" + str(start) + ".." + str(end))
		else:
			print(">" + seqid)
		seq=re.sub("(.{80})", "\\1\n", extractedSeq, 0, re.DOTALL)
		print(seq)


if outputPath != "":
	out.close()
