#!/usr/bin/env python
# -*- coding: utf8 -*-
#----------------------------------------------------------------------------
# Created By  : vloegler
# Created Date: 2022/01/04
# version ='1.0'
# ---------------------------------------------------------------------------
'''
This script retrieve contigs longer than X kb. 

It takes as input :
	-d --draft: a draft genome assembly to reorder (multi fasta)
	-m --minLength: minimum length of a contig (kb)
	-o --output: output prefix


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
# ---------------------------------------------------------------------------
# Definitions
def blastn(fasta1, fasta2, blastPath, out):
	blastCommand = blastPath + "blastn -query " + fasta1 + " -subject " + fasta2 + " -dust no -soft_masking false -max_target_seqs 1 -max_hsps 1 -outfmt '6 qseqid qlen length' -out " + out
	print(blastCommand)
	os.system(blastCommand)
# ---------------------------------------------------------------------------

# =============
# Get arguments
# =============

# Initiate the parser
parser = argparse.ArgumentParser()
parser.add_argument("-d", "--draft", help="draft genome assembly (multi fasta)", required=True)
parser.add_argument("-m", "--minLength", help="minimum length of a contig (kb)", type=int, default=1)
parser.add_argument("-o", "--output", help="Prefix of the output file", required=True)

# Read arguments from the command line
args = parser.parse_args()

draftPath=args.draft
minLength=args.minLength
outputPath=args.output

print("\n\t--- KEEPING NUCLEAR CONTIGS ---\n")
print("Arguments detected:")
print("\t--draft:\t"+draftPath)
print("\t--minLength:\t"+str(minLength))
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

# write output file
output=open(outputPath+".min" + str(minLength) + "kb.fasta", 'w')
kept = 0
for i in range(len(draftChr)):
	if draftLen[i] >= minLength * 1000:
		kept += 1
		seq=re.sub("(.{80})", "\\1\n", draftSeq[i], 0, re.DOTALL)
		output.write(draftLine[i])
		output.write(seq+"\n")
output.close()


print()
print(str(kept) + "/" + str(len(draftChr)) + " contigs kept")
print()
