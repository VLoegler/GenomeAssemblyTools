#!/usr/bin/env python
# -*- coding: utf8 -*-
#----------------------------------------------------------------------------
# Created By  : vloegler
# Created Date: 2022/09/15
# version ='1.4'
# ---------------------------------------------------------------------------
'''
This script extract a sequence from a fasta file

It takes as input :
	-f --fasta: fasta file
	-id --seqid: ID of the sequence, or multiple IDs separated by space
	-p --prefix: Prefix that the sequence ID have to match to be retained (several allowed seperated by space)
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
parser.add_argument("-id", "--seqid", help="ID of the sequence, or multiple IDs separated by space", nargs='+', default = [])
parser.add_argument("-p", "--prefix", help="Prefix that the sequence ID have to match to be retained (several allowed seperated by space)", nargs='+', default = [])
parser.add_argument("-s", "--start", help="Start position (inclusive, 1-based)", type = int, default = 0)
parser.add_argument("-e", "--end", help="End position (inclusive, 1-based)", type = int, default = 0)

# Read arguments from the command line
args = parser.parse_args()

fastaPath = args.fasta
outputPath=args.output
seqidList=args.seqid
prefixList=args.prefix
start=int(args.start)
end=int(args.end)

# Check if at least 1 seq ID or seq Prefix was given
if len(seqidList)+len(prefixList) == 0:
	raise ValueError("No sequence ID or Prefix provided. ")

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


# Check if all sequence ID or prefix are present in the fasta to extract
for seqid in seqidList:
	if seqid not in Chr:
		raise ValueError("Sequence ID " + seqid + " is not in fasta file.")
for prefix in prefixList:
	if not any(item.startswith(prefix) for item in Chr):
		raise ValueError("Prefix " + prefix + " is not in fasta file.")

for ID in Chr:
	# Check if ID is in the seqid list or starts with any prefix in prefixList
	if (ID in seqidList) or any(ID.startswith(item) for item in prefixList):

		# Retrieve sequence corresponding to the ID
		targetSeq = Seq[Chr.index(ID)]

		# If coordiantes are present, check proper coordinates 
		if start != 0 or end != 0:
			if start < 1 or start > len(targetSeq) or end < 1 or end > len(targetSeq) or end < start:
				raise ValueError("Bad coordinates. Start and End must be 1-based, not higher than sequence length and starts < end.")

		# If coordinates are present, extract the fragment
		if start != 0 or end != 0:
			extractedSeq = targetSeq[(start-1):end]
		else:
			extractedSeq = targetSeq

		# Write the sequence to the ouptupt (Stdout or file)
		if outputPath != "":
			if start != 0 or end != 0:
				out.write(">" + ID + "_" + str(start) + ".." + str(end) + "\n")
			else:
				out.write(">" + ID + "\n")
			seq=re.sub("(.{80})", "\\1\n", extractedSeq, 0, re.DOTALL)
			if seq.endswith('\n'):
				out.write(seq)
			else:
				out.write(seq+"\n")
		else:
			if start != 0 or end != 0:
				print(">" + ID + "_" + str(start) + ".." + str(end))
			else:
				print(">" + ID)
			seq=re.sub("(.{80})", "\\1\n", extractedSeq, 0, re.DOTALL)
			if seq.endswith('\n'):
				print(seq[:-1]) # Remove the last \n
			else:
				print(seq)


if outputPath != "":
	out.close()
