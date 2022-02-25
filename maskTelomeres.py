#!/usr/bin/env python
# -*- coding: utf8 -*-
#----------------------------------------------------------------------------
# Created By  : vloegler
# Created Date: 2022/02/24
# version ='1.0'
# ---------------------------------------------------------------------------
'''
This script masks telomeres (specified by a size) of a given assembly. 
Nucleotides are replaced by N. 

It takes as input :
	-f --fasta: fasta file
	-o --output: output file
	-k --kbToHide: Size of regions to hide at chromosome extremities(kb, default = 20)
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
parser.add_argument("-o", "--output", help="output file", required=True)
parser.add_argument("-k", "--kbToHide", help="Size of regions to hide at chromosome extremities(kb, default = 20)", type=int, default=20)

# Read arguments from the command line
args = parser.parse_args()

fastaPath = args.fasta
outputPath=args.output
kbToHide=args.kbToHide

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
		Chr += [line]
		if seq != "":
			Seq += [seq]
		seq=""
	else:
		seq += line.split("\n")[0]
Seq += [seq]
fasta.close()

hiddenSeq = []
for i in range(len(Seq)):
	hidden = "N"*kbToHide*1000 + Seq[i][kbToHide*1000:len(Seq[i])-(kbToHide*1000)] + "N"*kbToHide*1000
	hiddenSeq += [hidden]

out = open(outputPath, "w")
for i in range(len(Chr)):
	out.write(Chr[i])
	seq=re.sub("(.{80})", "\\1\n", hiddenSeq[i], 0, re.DOTALL)
	out.write(seq+"\n")
out.close()
