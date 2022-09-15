#!/usr/bin/env python
# -*- coding: utf8 -*-
#----------------------------------------------------------------------------
# Created By  : vloegler
# Created Date: 2022/01/04
# version ='1.0'
# ---------------------------------------------------------------------------
'''
This script retrieve nuclear contigs and discard small contigs. A nuclear 
contigs is defined as a contig that aligns on the reference nuclear genome
over more than X consecutive kb with blast. 

It takes as input :
	-d --draft: a draft genome assembly to reorder (multi fasta)
	-db --database: the blast database of nuclear chromosomes (done with makeblastdb)
	-r --reference: reference of the nuclear genome (multi fasta)
	-m --minLength: minimum length of a contig (kb)
	-a --minAlignedLength: minimum length of the alignment of the contig on the ref (kb)
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
# ---------------------------------------------------------------------------
# Definitions
def blastnDB(query, database, blastPath, out):
	blastCommand = blastPath + "blastn -query " + query + " -db " + database + " -dust no -soft_masking false -max_target_seqs 1 -max_hsps 1 -num_threads 4 -outfmt '6 qseqid qlen length' -out " + out
	print(blastCommand)
	os.system(blastCommand)
def blastnRef(query, reference, blastPath, out):
	blastCommand = blastPath + "blastn -query " + query + " -subject " + reference + " -dust no -soft_masking false -max_target_seqs 1 -max_hsps 1 -num_threads 4 -outfmt '6 qseqid qlen length' -out " + out
	print(blastCommand)
	os.system(blastCommand)
# ---------------------------------------------------------------------------

# =============
# Get arguments
# =============

# Initiate the parser
parser = argparse.ArgumentParser()
parser.add_argument("-d", "--draft", help="draft genome assembly (multi fasta)", required=True)
parser.add_argument("-db", "--database", help="the blast database of nuclear chromosomes (done with makeblastdb). Either reference or Database must be given", default = "")
parser.add_argument("-r", "--reference", help="referance of the nuclear genome (multifasta). Either reference or Database must be given", default = "")
parser.add_argument("-m", "--minLength", help="minimum length of a contig (kb)", type=int, default=1)
parser.add_argument("-a", "--minAlignedLength", help="minimum alignment length on the nuclear genome (kb)", type=int, default=2)
parser.add_argument("-o", "--output", help="Prefix of the output file", required=True)
parser.add_argument("-b", "--blastPath", help="Path to blast+ function, if not in path", type=str, default="")

# Read arguments from the command line
args = parser.parse_args()

draftPath=args.draft
dbPath=args.database
refPath=args.reference
minLength=args.minLength
minAlignedLength=args.minAlignedLength
outputPath=args.output
blast=args.blastPath
if blast != "" and not blast.endswith("/") :
	blast += "/"
if refPath == "" and dbPath == "":
	raise ValueError("You must either specify the path to a blastn database or to a reference in fasta format")
elif refPath != "" and dbPath != "":
	raise ValueError("You cannot specify both the path to a blastn database and to a reference in fasta format")



print("\n\t--- KEEPING NUCLEAR CONTIGS ---\n")
print("Arguments detected:")
print("\t--draft:\t"+draftPath)
if dbPath != "":
	print("\t--database:\t"+dbPath)
else:
	print("\t--reference:\t"+refPath)	
print("\t--minLength:\t"+str(minLength)+"kb")
print("\t--minAlignedLength:\t"+str(minAlignedLength)+"kb")
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


# Run Blastn of draft against ref
start = time.time()
if dbPath != "":
	blastnDB(draftPath, dbPath, blast, blastResultsPath)
else:
	blastnRef(draftPath, refPath, blast, blastResultsPath)	
end = time.time()
print("Alignment done: ran in "+str(round(end-start))+"s")

# Blast result analysis
# Keep contigs with length >= minLength and that aligned on nuclear ref genome
# for >= minAlignedLength consecutive kb
contigs2keep = []
blastResultsFile = open(blastResultsPath)
blastResults = csv.reader(blastResultsFile, delimiter="\t")
for row in blastResults:
	contig = row[0]
	contigLength = int(row[1])
	alignmentLength = int(row[2])
	if contigLength >= minLength * 1000 and alignmentLength >= minAlignedLength * 1000:
		contigs2keep += [contig]
blastResultsFile.close()
os.remove(blastResultsPath)


# write output file
output=open(outputPath+".Nuclear.fasta", 'w')
if len(contigs2keep) < len(draftChr) :
	output2=open(outputPath+".NonNuclear.fasta", 'w')
for i in range(len(draftChr)):
	if draftChr[i] in contigs2keep:
		seq=re.sub("(.{80})", "\\1\n", draftSeq[i], 0, re.DOTALL)
		output.write(draftLine[i])
		output.write(seq+"\n")
	else:
		seq=re.sub("(.{80})", "\\1\n", draftSeq[i], 0, re.DOTALL)
		output2.write(draftLine[i])
		output2.write(seq+"\n")		
output.close()

print()
print(str(len(contigs2keep)) + "/" + str(len(draftChr)) + " contigs kept as NUCLEAR CONTIGS")
print()
