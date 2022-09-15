#!/usr/bin/env python
# -*- coding: utf8 -*-
#----------------------------------------------------------------------------
# Created By  : vloegler
# Created Date: 2022/09/15
# version ='1.0'
# ---------------------------------------------------------------------------
'''
This script retrieve centromeres in an assembly and count the number of
contigs containing more than 1 centromere (probably falsly merged 
chromosomes). 

It takes as input :
	-d --draft: a draft genome assembly to reorder (multi fasta)
	-c --centromere: reference sequence of centromeres (multi fasta)
	-o --output: output prefix
	-b --blastPath: Path to blast+, if not in path

blast+ is used in this script. If blast+ is not in the path, the path can
be added to the variable blast. 

Output will be:
Assembly	NbMergedContigs
Draft1		0
Draft2		0
Draft3		1
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
def blastn(query, subject, blastPath, out):
	blastCommand = blastPath + "blastn -query " + query + " -subject " + subject + " -word_size 10 -outfmt 6 -qcov_hsp_perc 90 -dust no -soft_masking false -perc_identity 80 -out " + out
	os.system(blastCommand)
# ---------------------------------------------------------------------------

# =============
# Get arguments
# =============

# Initiate the parser
parser = argparse.ArgumentParser()
parser.add_argument("-d", "--draft", help="draft genome assemblies (multi fasta)", nargs='+', required=True)
parser.add_argument("-c", "--centromere", help="reference sequence of centromeres (multi fasta)", required = True)
parser.add_argument("-o", "--output", help="Prefix of the output file", default = "")
parser.add_argument("-b", "--blastPath", help="Path to blast+ function, if not in path", type=str, default="")

# Read arguments from the command line
args = parser.parse_args()

draftPaths=args.draft
nbDraft=len(draftPaths)
centroPath=args.centromere
outputPath=args.output
blast=args.blastPath
if blast != "" and not blast.endswith("/") :
	blast += "/"


# Write header in output file
if outputPath != "":
	out = open(outputPath, 'w')
	out.write("Assembly\tNbMergedContigs\n")
else:
	print("Assembly\tNbMergedContigs")


for d in range(nbDraft):
	draftName=draftPaths[d].split("/")[-1]
	draftPath = draftPaths[d]

	# Run blast of centromere against draft assembly
	blastResultsPath = datetime.now().strftime("%Y_%m_%d_%H_%M_%S_%m_%f") + "_" + str(randint(0, 10000)) + ".blastn"
	blastn(centroPath, draftPath, blast, blastResultsPath)	

	# Get all contigs
	draftChr = []
	draft=open(draftPath, 'r')
	for line in draft.readlines():
		if line.startswith(">"):
			draftChr += [line.strip().split(">")[1].split(" ")[0].split("\t")[0]]
	draft.close()

	# Count number of centromere per contig
	centromeres = [0]*len(draftChr)
	blastResults = open(blastResultsPath, "r")
	for line in blastResults.readlines():
		Chr = line.split("\t")[1]
		centromeres[draftChr.index(Chr)] += 1
	blastResults.close()
	os.remove(blastResultsPath)

	# Count number of contigs with more than 1 centromere
	count = 0
	for i in range(len(draftChr)):
		if centromeres[i] > 1:
			count += 1

	if outputPath != "":
		out.write(draftName + "\t" + str(count) + "\n")
	else:
		print(draftName + "\t" + str(count))

if outputPath != "":
	out.close()


