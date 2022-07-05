#!/usr/bin/env python
# -*- coding: utf8 -*-
#----------------------------------------------------------------------------
# Created By  : vloegler
# Created Date: 2022/02/26
# version ='1.0'
# ---------------------------------------------------------------------------
'''
This script contain all tools and classes to deal with BED coordinates and fasta
files. 
'''
# ---------------------------------------------------------------------------
#from __future__ import annotations
import os
import sys
import time
import csv
from typing import List
# ---------------------------------------------------------------------------

def rank_simple(vector):
	return sorted(range(len(vector)), key=vector.__getitem__)

def decreasing_rank_simple(vector):
	return sorted(range(len(vector)), key=vector.__getitem__)[::-1]
# ---------------------------------------------------------------------------
# Class

# First declare class names
class BEDcoordinates: pass
class BED: pass

class BEDcoordinates:
	def __init__(self, id:str, start:int, end:int):
		# Class constructor with reel coordinates
		if start < end:
			self.void = False
			self.id = id
			self.start = start
			self.end = end
		elif start == end:
			self.void = True
		else:
			raise Exception("Wrong coordinates given for BEDcoordinates object. End must be larger than start. ")
	def copy(self):
		return BEDcoordinates(self.id, self.start, self.end)
	def __str__(self):
		if self.void:
			return "Void"
		else:
			return self.id+"\t"+str(self.start)+"\t"+str(self.end)
	def overlap(self, b:BEDcoordinates):
		if self.void or b.void:
			overlap = False
		else:
			# Store b1+b2 in newBED object
			newBED=[]
			if self.id == b.id:
				if b.start <= self.start:
					if b.end <= self.start:
						overlap = False
					elif b.end < self.end:
						overlap = True
					else: 
						overlap = True
				elif b.start < self.end:
					if b.end < self.end:
						overlap = True
					else: 
						overlap = True
				else:
					overlap = False
			else: 
				overlap = False
		return overlap
	def overlapLen(self, b:BEDcoordinates):
		if self.void or b.void or self.id != b.id :
			overlap = 0
		elif b.start <= self.start:
			if b.end <= self.start:
				overlap = 0
			elif b.end < self.end:
				overlap = b.end - self.start
			else: 
				overlap = self.end - self.start
		elif b.start < self.end:
			if b.end < self.end:
				overlap = b.end - b.start
			else: 
				overlap = self.end - b.start
		else:
			overlap = 0
		return overlap
	def addCoordinates(self, b:BEDcoordinates):
		if self.void and b.void:
			newBEDcoordinates = []
		elif self.void:
			newBEDcoordinates = [b]
		elif b.void:
			newBEDcoordinates = [self]
		else:
			# Store b1+b2 in newb BED object
			if self.id == b.id:
				if b.start <= self.start:
					if b.end <= self.start:
						newBEDcoordinates = [b, self]
					elif b.end < self.end:
						newBEDcoordinates = [BEDcoordinates(id = self.id, start = b.start, end = self.end)]
					else: 
						newBEDcoordinates = [b] # b1 is comprised in b2
				elif b.start < self.end:
					if b.end < self.end:
						newBEDcoordinates = [self] # b2 is comprised in b1
					else: 
						newBEDcoordinates = [BEDcoordinates(id = self.id, start = self.start, end = b.end)]
				else:
					newBEDcoordinates = [self, b]
			else: 
				newBEDcoordinates = [self, b]
			# Return the addition as BED object (list of BED coordinates)
		return newBEDcoordinates
	def substractCoordinates(self, b:BEDcoordinates):
		if self.void:
			newBEDcoordinates = []
		elif b.void:
			newBEDcoordinates = [self]
		elif b.id == self.id:
			if b.start <= self.start:
				if b.end <= self.start:
					newBEDcoordinates = [self]
				elif b.end < self.end:
					newBEDcoordinates = [BEDcoordinates(id = self.id, start = b.end, end = self.end)]
				else: 
					newBEDcoordinates = [] # the BED coordinates are deleted
			elif b.start < self.end:
				if b.end < self.end:
					newBEDcoordinates = [BEDcoordinates(id = self.id, start = self.start, end = b.start), BEDcoordinates(id = self.id, start = b.end, end = self.end)]
				else: 
					newBEDcoordinates = [BEDcoordinates(id = self.id, start = self.start, end = b.start)]
			else:
				newBEDcoordinates = [self]
		else: 
			newBEDcoordinates = [self]
		# Return the substraction as BED object (list of BED coordinates)
		return newBEDcoordinates


class BED:
	def __init__(self, *args):
		if len(args) == 0:
			# Class constructor for void coordinates
			self.IDs = []
			self.coordinates = []
			self.len = 0
			self.nbIDs = 0
		else:
			# Class constructor for BED and BEDcoordinates arguments
			self.IDs = []
			self.coordinates = []
			coordinates2add = []
			for arg in args:
				if isinstance(arg, BEDcoordinates):
					if not arg.void:
						coordinates2add += [arg]
				elif isinstance(arg, BED):
					if not arg.nbIDs == 0:
						for i in range(len(arg.IDs)):
							coordinates2add += arg.coordinates[i]
				elif isinstance(arg, list):
					for x in arg:
						if isinstance(x, BEDcoordinates):
							if not x.void:
								coordinates2add += [x]
						elif isinstance(x, BED):
							if not x.nbIDs == 0:
								for i in range(len(x.IDs)):
									coordinates2add += x.coordinates[i]
						else:
							raise Exception("Wrong type argument given. Only BED and BEDcoordinates arguments taken. ")
				else:
					raise Exception("Wrong type argument given. Only BED and BEDcoordinates arguments taken. ")
				for b in coordinates2add:
					if b.id in self.IDs:
						self.coordinates[self.IDs.index(b.id)] += [b]
					else:
						self.IDs += [b.id]
						self.coordinates += [[b]]
			self.nbIDs = len(self.coordinates)
			if self.nbIDs == 0:
				self.IDs = []
				self.coordinates = []
				self.len = 0
				self.nbIDs = 0
			else:
				if max([len(x) for x in self.coordinates]) > 1: # If there is more than one BEDcoordinates per id
					self.removeOverlap()
				self.order()
				self.len = self.getLen()
	def removeOverlap(self):
		# This function check if there is an overlap between some coordinates
		# and merge them if it is the case
		for i in range(self.nbIDs):
			minPos = min([b.start for b in self.coordinates[i]])
			maxPos = max([b.end for b in self.coordinates[i]])
			totalBEDcoordinates = [BEDcoordinates(id = self.IDs[i], start = minPos, end = maxPos)]
			for b2substract in self.coordinates[i]:
				newTotalBEDcoordinates = []
				for b in totalBEDcoordinates:
					newTotalBEDcoordinates += b.substractCoordinates(b2substract)
				totalBEDcoordinates = newTotalBEDcoordinates.copy()
			# invert the obtained BED
			# Get all positions and order them
			allPositions = [minPos, maxPos]
			for b in totalBEDcoordinates:
				allPositions += [b.start, b.end]
			allPositions.sort()
			# 
			finalBEDcoordinates = []
			for j in range(0, len(allPositions), 2):
				finalBEDcoordinates += [BEDcoordinates(id = self.IDs[i], start = allPositions[j], end = allPositions[j+1])]
			self.coordinates[i] = finalBEDcoordinates.copy()
	def order(self):
		# This function order the BED by id then by start position
		# First order IDs
		alphabeticalIDs = sorted(self.IDs, key = str)
		indices = [self.IDs.index(id) for id in alphabeticalIDs]
		newCoordinates = [self.coordinates[i] for i in indices]
		self.IDs = alphabeticalIDs.copy()
		self.coordinates = newCoordinates.copy()
		# Secondly, order each BEDcoordinates
		for i in range(self.nbIDs):
			startPositions = [b.start for b in self.coordinates[i]]
			subOrder = rank_simple(startPositions)
			newCoordinates = [self.coordinates[i][j] for j in subOrder]
			self.coordinates[i] = newCoordinates.copy()
	def getLen(self):
		lengths = []
		for i in range(self.nbIDs):
			lengths += [b.end-b.start for b in self.coordinates[i]]
		return sum(lengths)
	def copy(self):
		return BED(self)
	def __str__(self):
		if self.nbIDs == 0:
			return "Void"
		else:
			toPrint=""
			for i in range(self.nbIDs):
				for b in self.coordinates[i]:
					toPrint +=  "[" + b.__str__() + "]\n"
			return toPrint[0:-1]
	def __add__(self, B):
		if isinstance(B, BED):
			return BED(self, B)
		else: 
			raise Exception("Wrong type argument given. Add operator only takes BED objects. ")
	def __sub__(self, B):
		if isinstance(B, BED):
			newBED = self.copy()
			newBED.substractBED(B)
			return newBED
		else: 
			raise Exception("Wrong type argument given. Sub operator only takes BED objects. ")
	#def getCoordinates(self, index:int):
	#	return self.coordinates[index]
	def getID(self, id:str):
		# return a BED of all coordinates with id
		if id in self.IDs:
			return BED(self.coordinates[self.IDs.index(id)])
		else:
			return BED()
	#def removeCoordinates(self, index:int):
	#	newCoordinates = self.coordinates.copy()
	#	del newCoordinates[index]
	#	self.coordinates = newCoordinates
	#	self.nbIDs = len(self.coordinates)
	#	self.len = self.getLen()
	def addCoordinates(self, b:BEDcoordinates):
		if b.void:
			pass
		elif self.nbIDs == 0:
			self.IDs = [b.id]
			self.coordinates = [[b]]
		else:
			if b.id not in self.IDs:
				self.IDs += [b.id]
				self.coordinates += [[b]]
			else:
				i = self.IDs.index(b.id)
				overlap = False
				overlappedCoordinates = []
				coordinates2del = []
				for j in range(len(self.coordinates[i])):
					if self.coordinates[i][j].overlap(b):
						overlap = True
						overlappedCoordinates += [self.coordinates[i][j]]
						coordinates2del += [j]
				if not overlap:
					self.coordinates[i] +=  [b]
				else:
					newCoordinates = self.coordinates[i].copy()
					# Remove Overlapped coordinates from BED
					coordinates2del.sort(reverse = True) # Sort reverse to remove last indexes first
					for j in coordinates2del:
						del newCoordinates[j]
					# Sequentially add Overlapped coordinates to b
					for j in range(len(overlappedCoordinates)):
						b = b.addCoordinates(overlappedCoordinates[j])[0]
					self.coordinates[i] = newCoordinates + [b]
					self.order()
		self.nbIDs = len(self.coordinates)
		self.len = self.getLen()
	def checkOverlap(self):
		# This function check if there is an overlap between coordinates and return True of False
		for i in range(self.nbIDs):
			for j1 in range(len(self.coordinates[i])):
				for j2 in range(len(self.coordinates[i])):
					if j1 < j2:
						if self.coordinates[i][j1].overlap(self.coordinates[i][j2]):
							return True
		return False
	def substractCoordinates(self, b2:BEDcoordinates):
		if self.nbIDs == 0 or b2.void:
			pass
		else:
			if b2.id not in self.IDs:
				pass
			else: 
				# Remove b2 to each coordinates with same ID
				indexId =  self.IDs.index(b2.id)
				newCoordinates = []
				for b1 in self.coordinates[indexId]:
					newCoordinates += b1.substractCoordinates(b2)
				self.coordinates[indexId] = newCoordinates
				self.len = self.getLen()
	def substractBED(self, B:BED):
		if B.nbIDs == 0:
			pass
		else:
			# Remove each coordinates of BED to this BED object
			for i in range(len(B.IDs)):
				for b in B.coordinates[i]:
					self.substractCoordinates(b)
	def overlapLen(self, B:BED, percent = False):
		overlap = 0
		for id in B.IDs:
			if id in self.IDs:
				for b1 in self.coordinates[self.IDs.index(id)]:
					for b2 in B.coordinates[B.IDs.index(id)]:
						overlap += b1.overlapLen(b2)
		if percent :
			if self.len == 0:
				raise Exception("Cannot compute overlap length percent on a BED with length 0.")
			else:
				return (overlap / self.len) * 100
		else :
			return overlap


# ---------------------------------------------------------------------------
# Definitions
def getNonNBED(seqName:str, seq:str):
	indices = [index for index, element in enumerate(seq) if (element == "n" or element == "N")]
	# Get BED positions of Ns:
	coordinates = []
	for i in range(len(indices)):
		if i == 0: # If first N
			if indices[i+1] == indices[i]+1:
				coordinates += [indices[i]] # Add starting coordinates
			else: # First N in the sequence is a single N
				coordinates += [indices[i]] # Add starting coordinate
				coordinates += [indices[i]+1] # Add ending coordinate
		elif i == len(indices)-1: # if last N
			if indices[i-1] == indices[i]-1:
				coordinates += [indices[i]+1] # Add ending coordinate
			else : # Last N in the sequence is a single N
				coordinates += [indices[i]] # Add starting coordinate
				coordinates += [indices[i]+1] # Add ending coordinate
		else:
			cond1 = (indices[i-1] == indices[i]-1)
			cond2 = (indices[i+1] == indices[i]+1)
			if not cond1 and not cond2: # Case of a single n in the sequence
				coordinates += [indices[i]] # Add starting coordinate
				coordinates += [indices[i]+1] # Add ending coordinate
			elif (cond1 ^ cond2): # XOR operator
				if indices[i-1] == indices[i]-1:
					coordinates += [indices[i]+1] # Add ending coordinate
				else :
					coordinates += [indices[i]] # Add starting coordinate
	N_BED = BED()
	for i in range(0, len(coordinates), 2):
		startPos = coordinates[i] + 1
		endPos = coordinates[i+1] + 1
		N_BED.addCoordinates(BEDcoordinates(id = seqName, start = startPos, end = endPos))
	# Invert Bed Coordinates
	# First create BED with all position
	NonN_BED = BED(BEDcoordinates(id = seqName, start = 1, end = len(seq)+1))
	NonN_BED.substractBED(N_BED) # Substract N Positions
	return NonN_BED

def getBED(fastaPath):
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
	# Create the BED file for each chromosome [ChrID, Start, End]
	# Start is included, End is excluded
	BED2Add = []
	fastaBED = BED()
	for i in range(len(Chr)):
		# Get non N positions
		BED2Add += [getNonNBED(Chr[i], Seq[i])]
	# Merging all BEDs into One
	fastaBED = BED(BED2Add)
	return fastaBED


'''
def rank_simple(vector):
	return sorted(range(len(vector)), key=vector.__getitem__)

def decreasing_rank_simple(vector):
	return sorted(range(len(vector)), key=vector.__getitem__)[::-1]

def getShowCoords(refPath, draftPath, prefix, mummerPath, threads):
	# Run nucmer
	#nucmerShellCommand=mummerPath+"nucmer -t "+str(threads)+" --maxmatch --prefix "+prefix+" "+refPath+" "+draftPath
	nucmerShellCommand=mummerPath+"nucmer -t "+str(threads)+" --mum --prefix "+prefix+" "+refPath+" "+draftPath
	os.system(nucmerShellCommand)
	# Run show coords
	showcoordsShellCommand=mummerPath+"show-coords -TH "+prefix+".delta > "+prefix+".coords"
	os.system(showcoordsShellCommand)
	# Remove nucmer delta file
	os.remove(prefix+".delta")

'''
# ---------------------------------------------------------------------------

if __name__ == "__main__":
	print("Empty bed")
	a = BED()
	print(a)
	print("Simple non overlapping BED")
	b = BED(BEDcoordinates("a", 2, 100), BEDcoordinates("a", 300, 600))
	print(b)
	print("Overlapping bed")
	c = BED(BEDcoordinates("b", 2, 400), BEDcoordinates("b", 300, 600))
	print(c)
	print("Add coordinates")
	c.addCoordinates(BEDcoordinates("a", 500, 610))
	print(c)
	print("Merge 3 BEDs")
	d = BED(a, b, c)
	print(d)
	print("getLen()")
	print(d.getLen())
	print("Create BED with BED, BEDcoordinates and list of BEDcoordinates objects")
	e = BED(c, b, BEDcoordinates(id = "aa", start = 3, end = 7), BEDcoordinates(id = "z", start = 3, end = 7), [BEDcoordinates(id = "z", start = 30, end = 70), BEDcoordinates(id = "z", start = 5, end = 15)])
	print(e)
	print("substractCoordinates")
	e.substractCoordinates(BEDcoordinates("aa", start = 4, end = 6))
	print(e)
	print("substractBED")
	e.substractBED(BED(BEDcoordinates(id = "aa", start = 1, end = 543), BEDcoordinates(id = "z", start = 1, end = 40)))
	print(e)
	print("Add operator")
	print("a")
	a = BED(BEDcoordinates("a", 2, 300))
	print(a)
	print("b")
	b = BED(BEDcoordinates("a", 200, 400))
	print(b)
	print("a+b")
	c = a+b
	print(c)
	print("Sub operator")
	print("a")
	a = BED(BEDcoordinates("a", 2, 100))
	print(a)
	print("b")
	b = BED(BEDcoordinates("a", 20, 41))
	print(b)
	print("a-b")
	c = a-b
	print(c)

	# Create huge BED
	import random
	import time
	id = "HUGE"
	positions = list(range(1, 10001))
	random.Random(1).shuffle(positions)
	coordinates = []
	for i in range(0, len(positions), 2):
		if positions[i] < positions[i+1]:
			startPos = positions[i]
			endPos = positions[i+1]
		else:
			startPos = positions[i+1]
			endPos = positions[i]
		coordinates += [BEDcoordinates(id = id, start = startPos, end = endPos)]
	print("Start merging HUGE BED")
	startTime = time.time()
	a = BED(coordinates)
	endTime = time.time()
	print(a)
	print("Ran in " + str(round(endTime - startTime)) + "s")


