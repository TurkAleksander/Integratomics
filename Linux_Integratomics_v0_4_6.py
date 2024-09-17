# TODOs
# Make the database use a new sqlite instance for each run of the analyses
# Import at least some additional feature to position conversion methods
# Implement option to disable nested permutations


# Import CherryPy global namespace
#import cherrypy 
#import web
import os
#import Cfg
import sqlite3 as lite
import sys
import numpy as np
import multiprocessing as mp
import scipy as sp
import os 
import glob
from random import choice
from matplotlib import pyplot as plt
from scipy.stats import gaussian_kde
import re
import time
import json
import shutil
from traceback import print_stack
import math
import shutil


#Class definitions
#This class takes vector with numbers representing numbers of genes in a region and returns the vector of blocks, of the same length as the input vector. This function requires minBlockProp value as an input, to prevent excessive fragmentation of the blocks. 
class makeBlocks:
	def __init__(self, genes, minBlockProp):
		self.genes = genes
		self.objects = sorted(unique(self.genes))
		self.minBlockProp = minBlockProp
		self.blocks = self.getBlocks()
		
	def getBlocks(self):
		for index, obj in enumerate(self.objects):
			if (len(self.genes[self.genes == obj])/float(len(self.genes)) < self.minBlockProp) and (1+index<len(self.objects)):
				self.genes[self.genes == obj] = self.objects[index + 1]
			elif (len(self.genes[self.genes == obj])/float(len(self.genes)) < self.minBlockProp) and (1+index==len(self.objects)):
				self.genes[self.genes == obj] = self.objects[index - 1]
			else:
				next
		return self.genes

#Main permutation class
class perMat:
		def __init__(self,mat,studyW,perm,genes):
			
			self.mat = mat
			self.nrow = np.shape(self.mat)[0]
			if np.ndim(self.mat)==1:
				self.mat = np.reshape(self.mat, (self.nrow, 1))
			self.ncol = np.shape(self.mat)[1]
			self.perm = perm
			self.studyW = studyW
			self.genes = genes
			self.blocks = makeBlocks(self.genes, 0.1).blocks
			self.FDRs = []
			self.FD = np.zeros(self.nrow)
			self.makeRR()

		def makeRR(self):
			self.matRR = np.zeros((self.ncol,self.nrow)).T
			for i in range(0,self.ncol):
				self.matRR[:,i]=((np.array(rankdata_rand(self.mat[:,i])))/float(self.nrow))**self.studyW[i]

			self.dataProd = np.prod(self.matRR, axis=1)
			self.dataRank = np.array(rankdata(self.dataProd))
			return(self.matRR, self.dataRank, self.dataProd)
		
		def permutationTest(self):
			
			#for i in range(0,self.perm):
			self.matRRperm = self.matRR
			for j in range(0,self.ncol):
				for block in unique(self.blocks):
					tempvec = self.matRRperm[self.blocks==block,j]
					np.random.shuffle(tempvec)
					self.matRRperm[self.blocks==block,j] = tempvec
					#np.random.shuffle(self.matRR[:,j]) #This will omit block like shuffling
			temp = np.array((self.dataProd,np.prod(self.matRRperm, axis=1)))
			temp = np.reshape(temp, 2*self.nrow)
			tempRanks = rankdata_max(temp)
			tempRanks=np.reshape(tempRanks,(2,self.nrow)).T
			self.FD += np.reshape(tempRanks[:,0]-self.dataRank,self.nrow)
			#print np.min(self.FD)

	#		self.FDRs = FD/float(self.perm)
	#		return(self.FDRs)


		def getScores(self):
			if self.ncol > 1:
				print("Many substudies!")
				try:
					self.FD[self.FD==0] = 1
					self.FDRs = self.FD/float(self.nrow*self.perm)
				except:
					print("Something is wrong!")
				self.scores = -np.log10(self.FDRs)
			elif self.ncol == 1:
				print ("Only one substudy, not doing permutations, passing scores to final prioritization...")
				self.scores = np.reshape(self.mat, self.nrow)
			else:
				print ("Warning: Problems with input matrix dimensions!")
			#return self.scores
			
		def getSimpleScores(self):
			if self.ncol > 1:
				print ("Calculating rank product for multiple substudies...")
				self.simplescores = np.prod(self.matRR, axis=1)
				self.simplescores = -np.log10(self.simplescores)
			elif self.ncol == 1:
				print ("Only one substudy, passing directly to second tier...")
				self.simplescores = np.reshape(self.mat, self.nrow)
			else:
				print ("Warning: Problems with input matrix dimensions!")
			#return self.simplescores

print ("Ales")

class GWplot:
	def __init__(self, idxs, chrs, scores, sid):
		self.chrs = chrs
		self.scores = scores
		self.idxs = idxs
		self.sid = sid
		
		self.chromosomes = np.array(unique_sorted(self.chrs))
		self.ticks=np.zeros(len(self.chromosomes))

		plt.figure(figsize=(15,6.5))
		
		for i in np.arange(len(self.chromosomes)):
			if i % 2:
				#print self.chromosomes[i]
				self.ticks[i] = np.mean(self.idxs[self.chrs == self.chromosomes[i]])
				plt.scatter(self.idxs[self.chrs==self.chromosomes[i]], self.scores[self.chrs==self.chromosomes[i]], c="r")
			else:
				#print "Other %s" % self.chromosomes[i] 
				#print self.idxs[self.chrs == self.chromosomes[i]]
				self.ticks[i] = np.mean(self.idxs[self.chrs == self.chromosomes[i]])
				plt.scatter(self.idxs[self.chrs==self.chromosomes[i]], self.scores[self.chrs==self.chromosomes[i]], c="b")
		print (self.ticks)
		print (self.chromosomes)
		plt.xticks(self.ticks, self.chromosomes, size='small', rotation=30)
		plt.xlabel("Chromosomal position")
		plt.ylabel("Integratomic scores")
		plt.ax = plt.gca()
		plt.ax.spines['right'].set_color('none')
		plt.ax.spines['top'].set_color('none')
		plt.ax.xaxis.set_ticks_position('bottom')
		plt.ax.yaxis.set_ticks_position('left')
		plt.ax.autoscale(enable=True, axis='both', tight=True)
		plt.ylim((-0.1,np.max(self.scores)+0.5))

		self.figfile = 'C:/Users/gkgen69/Desktop/Aleksander/Integratomics MS/Ales_data/MS_Integratomika_2024/Koda/Integratomika/integr/website/static/uploadify/wgplots/%s.png' % self.sid
		plt.savefig(self.figfile, bbox_inches=0)
		plt.clf()
		plt.close()
		
#Function definitions
def precInt(start, interval, overlap):
	start = int(start)
	jump = interval-overlap
	intStart = start - (start % jump)
	return intStart
	
def flatten(lst):
    result = []
    for element in lst: 
        if hasattr(element, '__iter__'):
            result.extend(flatten(element))
        else:
            result.append(element)
    return result

def rank_simple(vector):
    return sorted(range(len(vector)), key=vector.__getitem__)

def rankdata_dec(a):
	n = len(a)
	ivec = a.argsort()
	svec = np.sort(a)
	sumranks = 0
	dupcount = 0
	newarray = [0]*n
	for i in xrange(n):
		sumranks += i
		dupcount += 1
		if i==n-1 or svec[i] != svec[i+1]:
			botrank = n-i+dupcount-1
			for j in xrange(i-dupcount+1,i+1):
				newarray[ivec[j]] = botrank
			sumranks = 0
			dupcount = 0
	return newarray


def rankdata(a):
	n = len(a)
	ivec = a.argsort()
	svec = np.sort(a)
	sumranks = 0
	dupcount = 0
	newarray = [0]*n
	for i in xrange(n):
		sumranks += i
		dupcount += 1
		if i==n-1 or svec[i] != svec[i+1]:
			botrank = i + 1
			for j in xrange(i-dupcount+1,i+1):
				newarray[ivec[j]] = botrank
			sumranks = 0
			dupcount = 0
	return newarray

def rankdata_max(a):
	n = len(a)
	ivec = a.argsort()
	svec = np.sort(a)
	sumranks = 0
	dupcount = 0
	newarray = [0]*n
	for i in xrange(n):
		sumranks += i
		dupcount += 1
		if i==n-1 or svec[i] != svec[i+1]:
			botrank = i - dupcount + 2
			for j in xrange(i-dupcount+1,i+1):
				newarray[ivec[j]] = botrank
			sumranks = 0
			dupcount = 0
	return newarray

def rankdata_test(a):
	n = len(a)
	ivec = a.argsort()
	svec = np.sort(a)
	sumranks = 0
	dupcount = 0
	newarray = [0]*n
	for i in xrange(n):
		sumranks += i
		dupcount += 1
		if i==n-1 or svec[i] != svec[i+1]:
			botrank = np.arange((i-dupcount+2),(i+2),1)
			for j in xrange(i-dupcount+1,i+1):
				newarray[ivec[j]] = choice(botrank)
			sumranks = 0
			dupcount = 0
	return newarray

def clearZeros(vec):
	for i in range(len(vec)):
		if (vec[i] == None):
			vec[i] = 0
			return vec
	
def unique(seq):
   # Not order preserving
   keys = {}
   for e in seq:
       keys[e] = 1
   return keys.keys()

def unique_sorted(seq): 
   # order preserving
   checked = []
   for e in seq:
       if e not in checked:
           checked.append(e)
   return checked

def rankdata_rand(a):
	n = len(a)
	ivec = a.argsort()
	svec = np.sort(a)
	sumranks = 0
	dupcount = 0
	newarray = [0]*n
	for i in xrange(n):
		sumranks += i
		dupcount += 1
		if i==n-1 or svec[i] != svec[i+1]:
			botrank = np.arange((n-i+dupcount-1),(n-i-1),-1)
			for j in xrange(i-dupcount+1,i+1):
				newarray[ivec[j]] = choice(botrank)
			sumranks = 0
			dupcount = 0
	return newarray



def ales():
	yield "Permutation  <br/>"
	for i in range(10):
		yield "Permutation %s <br/>" % i
		#print i
		
		
#################################################################
######## IMPORT STUDY TYPE FROM STUDIES
#################################################################

class getStudyTypes:
	def __init__(self, sid):
		self.studyFiles = []
		self.studyTypes = []
		self.studyNames = []
		self.annoTypes=[]
		self.sid = sid

		# Get weights for study types
		studyloc = 'C:/Users/gkgen69/Desktop/Aleksander/Integratomics MS/Ales_data/MS_Integratomika_2024/Koda/Integratomika/integr/temp/%s/Study*.txt' % sid
		for name in glob.glob(studyloc):
			z=open(name, 'r')
			self.studyTypes.append(re.sub(r'\W+', '', z.readline().rstrip().replace('-', '').replace(' ', '')))
			self.studyNames.append(re.sub(r'\W+', '', z.readline().rstrip().replace('-', '').replace(' ', '')))
			self.annoTypes.append(re.sub(r'\W+', '', z.readline().rstrip()))
			z.close()	
			
		# Detect syudyTypes
		self.uniqStudyTypes = unique(self.studyTypes)
		self.typeW = []
		
def unique(seq):
   # Not order preserving
   keys = {}
   for e in seq:
       keys[e] = 1
   return keys.keys()
		
		
def importAnnot():
	annotFiles = []
	annotTypes = []
	annotNames = []

	annoloc = 'C:/Users/gkgen69/Desktop/Aleksander/Integratomics MS/Ales_data/MS_Integratomika_2024/Koda/Integratomika/integr/annot/*.txt'
	
	for annot in glob.glob(annoloc):

		z=open(annot, 'r')

		annoType = []
		z.readline()
		annoType = z.readline()
		header = z.readline().rstrip().split('\t')

		print ("This annotation is %s" % annoType)
		all=[]

		for line in z.readlines():
			cols = line.rstrip().split('\t')
			if len(cols) == 4:
				if (cols[header.index("name")]=="" or cols[header.index("chr")]==""):
					next
				else:
					if re.match("^.*chr.*$", cols[header.index("chr")]):
						all.append((cols[header.index("name")], cols[header.index("chr")], cols[header.index("start")], cols[header.index("stop")]))
					else:
						all.append((cols[header.index("name")], "chr"+cols[header.index("chr")], cols[header.index("start")], cols[header.index("stop")]))
			else:
				next
		z.close()

		con = lite.connect('C:/Users/gkgen69/Desktop/Aleksander/Integratomics MS/Ales_data/MS_Integratomika_2024/Koda/Integratomika/integr/sqlite3/annotation.db')
		cur = con.cursor()

		sql0 = "DROP TABLE IF EXISTS %s" % annoType 
		cur.execute(sql0)
		sql1 = "CREATE TABLE %s (name TEXT, chr TEXT, start INT, stop INT)" % annoType 
		cur.execute(sql1)
		sql2 = "INSERT INTO %s VALUES(?, ?, ?, ?)" % annoType
		cur.executemany(sql2, all)
		con.commit()


def prepareBackbone(interval, overlap):
	z=open('C:/Users/gkgen69/Desktop/Aleksander/Integratomics MS/Ales_data/MS_Integratomika_2024/Koda/Integratomika/integr/data/genome_coor_hg19.txt', 'r')

	chrNames=[]
	chrLengths=[]
	studyType = "Backbone"

	for line in z.readlines():
		cols = line.split('\t')
		chrNames.append(cols[0])
		chrLengths.append(cols[1])
	z.close()

	all=[]
	idx = 0
	for j in range(len(chrLengths)):
		print ("Chromosome length %s <br />" % chrLengths[j])

		i = 1	
		chrLength = int(chrLengths[j])
		chrName = chrNames[j]
		i=interval-overlap
		while i + interval + overlap < chrLength:
			i = i + interval - overlap
			all.append((idx, chrName, i-interval, i-1))
			idx = idx + 1

	numIntervals = idx
	con = lite.connect('C:/Users/gkgen69/Desktop/Aleksander/Integratomics MS/Ales_data/MS_Integratomika_2024/Koda/Integratomika/integr/sqlite3/test.db')
	cur = con.cursor()

	sql0 = "DROP TABLE IF EXISTS %s" % studyType 
	cur.execute(sql0)
	sql1 = "CREATE TABLE %s (idx INT, chr TEXT, start INT, stop INT)" % studyType 
	cur.execute(sql1)
	sql2 = "INSERT INTO %s VALUES(?, ?, ?, ?)" % studyType
	cur.executemany(sql2, all)
	
	con.commit()

	numRegions = idx
	return(numRegions)
	


def hgncImport(interval, overlap, numRegions):
	z=open('C:/Users/gkgen69/Desktop/Aleksander/Integratomics MS/Ales_data/MS_Integratomika_2024/Koda/Integratomika/integr/data/Genes.txt', 'r')

	annoType = []
	z.readline()
	annoType = z.readline()
	header = z.readline().rstrip().split('\t')

	print ("This annotation is %s" % annoType)
	all=[]
	n=0

	for line in z.readlines():
		cols = line.rstrip().split('\t')
		if len(cols) == 4:
			if (cols[header.index("name")]=="" or cols[header.index("chr")]==""):
				next
			else:
				if not re.match("^.*chr.*$", cols[header.index("chr")]):
					cols[header.index("chr")] = "chr"+cols[header.index("chr")]
				if cols[header.index("start")]>cols[header.index("stop")]:
					cols[header.index("stop")], cols[header.index("start")] = cols[header.index("start")], cols[header.index("stop")]
				for i in range(precInt(cols[header.index("start")], interval, overlap), precInt(cols[header.index("stop")], interval, overlap)+1, interval-overlap):
					all.append((cols[header.index("chr")], i, 1, cols[header.index("name")]))
					n += 1
		else:
			next


	z.close()
	#yield "%s entries imported..." % idx

	con = lite.connect('C:/Users/gkgen69/Desktop/Aleksander/Integratomics MS/Ales_data/MS_Integratomika_2024/Koda/Integratomika/integr/sqlite3/test.db')
	cur = con.cursor()

	sql0 = "DROP TABLE IF EXISTS %s" % "Genes"
	cur.execute(sql0)
	sql1 = "CREATE TABLE %s (chr TEXT, start INT, score INT, name TEXT)" % "Genes"
	cur.execute(sql1)
	sql2 = "INSERT INTO %s VALUES(?, ?, ?, ?)" % "Genes"
	cur.executemany(sql2, all)

	con.commit()
	sql3 = "SELECT sum(%s.score) FROM Backbone LEFT OUTER JOIN %s ON %s.start = Backbone.start AND Backbone.chr = %s.chr GROUP BY Backbone.idx" % ("Genes", "Genes", "Genes", "Genes")
	cur.execute(sql3)
	output = cur.fetchall()
	genes = clearZeros(np.reshape(output, numRegions))

	return(genes)


def importCoorStudy(studyFiles, study, studyTypes, studyNames, interval, overlap, cur, con):
	print ("This are coordinates information for %s:" % studyFiles[study])
	z=open(studyFiles[study], 'r')

	z.readline()
	z.readline()
	z.readline()

	print ("This study is on %s, done by %s" % (studyTypes[study], studyNames[study]))
	
	chr=[]
	start=[]
	stop=[]
	score=[]
	all=[]
	output=[]
	NAval = 0

	for line in z.readlines():
		cols = re.split('\t| |\n', line.rstrip())
		cols[3]=float(cols[3].replace(',','.'))
		try:
			int(cols[1])
			if cols[1]>cols[2]:
				for i in range(precInt(cols[2], interval, overlap), precInt(cols[1], interval, overlap)+1, interval-overlap):
					all.append((cols[0], i, cols[3]))
				#print "Start position %s was downstream of stop position %s on chromosome %s! Using stop position as start..." %  (cols[1], cols[2], cols[0])
			else:
				for i in range(precInt(cols[1], interval, overlap), precInt(cols[2], interval, overlap)+1, interval-overlap):
					all.append((cols[0], i, cols[3]))
				score.append(cols[3])
		except:
			NAval += 1
	    #all.append((cols[0:4]))
	#yield 'Removed %s missing values<br/>' % NAval
	z.close()

	if (max([x[2] for x in all])<1):
		print ("This does not seem to be a dataset in Scores, but P values - so converting by -log10 transformation...")
		for i,e in enumerate(all):
			temp=list(all[i])
			print (temp[2])
			temp[2]=-math.log10(temp[2])
			all[i]=tuple(temp)
	else:
		next

	#print idx, "entries imported..."

	sql0 = "DROP TABLE IF EXISTS %s" % studyNames[study]
	cur.execute(sql0)
	sql1 = "CREATE TABLE %s (chr TEXT, start INT, score INT)" % studyNames[study]
	cur.execute(sql1)
	sql2 = "INSERT INTO %s VALUES(?, ?, ?)" % studyNames[study]
	cur.executemany(sql2, all)

	con.commit()
	sql3 = "SELECT avg(%s.score) FROM Backbone LEFT OUTER JOIN %s ON %s.start = Backbone.start AND Backbone.chr = %s.chr GROUP BY Backbone.idx" % (studyNames[study], studyNames[study], studyNames[study], studyNames[study])
	cur.execute(sql3)
	output = cur.fetchall()

	return output


def importAnnotStudy(studyFiles, study, studyTypes, studyNames, annoTypes, interval, overlap, cur, con):
	#print "The output is in the form of gene names"
	all=[]

	z=open(studyFiles[study], 'r')

	z.readline()
	z.readline()
	z.readline()

	for line in z.readlines():
		cols = re.split('\t| |\n', line.rstrip())
		if len(cols[0])>0:
			all.append((cols[0], float(cols[1].replace(',','.'))))
		else:
			continue
	z.close()
	print (all)
	
	if (max([x[1] for x in all])<1):
		print ("This does not seem to be a dataset in Scores, but P values - so converting by -log10 transformation...")
		for i,e in enumerate(all):
			temp=list(all[i])
			print (temp[1])
			temp[1]=-math.log10(temp[1])
			all[i]=tuple(temp)
	else:
		next

	con = lite.connect('C:/Users/gkgen69/Desktop/Aleksander/Integratomics MS/Ales_data/MS_Integratomika_2024/Koda/Integratomika/integr/sqlite3/test.db')
	cur = con.cursor()

	sql0 = "DROP TABLE IF EXISTS %s" % studyNames[study] 
	print (sql0)
	cur.execute(sql0)
	sql1 = "CREATE TABLE %s (name TEXT, score INT)" % studyNames[study] 
	cur.execute(sql1)
	sql2 = "INSERT INTO %s VALUES(?, ?)" % studyNames[study]
	cur.executemany(sql2, all)
	con.commit()

	cur.execute("ATTACH 'C:/Users/gkgen69/Desktop/Aleksander/Integratomics MS/Ales_data/MS_Integratomika_2024/Koda/Integratomika/integr/sqlite3' AS annotation")
	sql3 = "SELECT chr, start, stop, score FROM annotation.%s INNER JOIN %s ON annotation.%s.name = %s.name" % (annoTypes[study], studyNames[study], annoTypes[study], studyNames[study])
	cur.execute(sql3)
	output = cur.fetchall()
	all = []
	for cols in output:
		if cols[1]>cols[2]:
			for i in range(precInt(cols[2], interval, overlap), precInt(cols[1], interval, overlap)+1, interval-overlap):
				all.append((cols[0], i, cols[3]))
			#print "Start position %s was downstream of stop position %s on chromosome %s! Using stop position as start..." %  (cols[1], cols[2], cols[0])
		else:
			for i in range(precInt(cols[1], interval, overlap), precInt(cols[2], interval, overlap)+1, interval-overlap):
				all.append((cols[0], i, cols[3]))
		#score.append(cols[3])


	sql0 = "DROP TABLE IF EXISTS %s" % studyNames[study]
	cur.execute(sql0)
	sql1 = "CREATE TABLE %s (chr TEXT, start INT, score INT)" % studyNames[study]
	cur.execute(sql1)
	sql2 = "INSERT INTO %s VALUES(?, ?, ?)" % studyNames[study]
	cur.executemany(sql2, all)


	con.commit()
	sql3 = "SELECT avg(%s.score) FROM Backbone LEFT OUTER JOIN %s ON %s.start = Backbone.start AND Backbone.chr = %s.chr GROUP BY Backbone.idx" % (studyNames[study], studyNames[study], studyNames[study], studyNames[study])
	#sql3 = "SELECT sum(%s.score) FROM Backbone LEFT OUTER JOIN %s ON %s.start = Backbone.start AND Backbone.chr = %s.chr GROUP BY Backbone.idx" % (studyNames[study], studyNames[study], studyNames[study], studyNames[study])
	cur.execute(sql3)
	output = cur.fetchall()
	
	return output

	print ("Importing study %s complete" % studyNames[study])

