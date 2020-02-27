import sys
import csv
import os
import re
import numpy as np
import argparse
import statsmodels
from statsmodels.stats import multitest
from scipy.stats import ttest_ind
from scipy.stats import variation
from scipy.stats.mstats import gmean
from scipy.stats import rankdata
import pandas as pd
import random
import math
from tkinter import *
import tkinter as tk



def readFile(infile):
	global firstRow
	fn = infile
	with open(fn, 'r') as myfile:
	    file = csv.reader(myfile, dialect='excel')
	    
	    d = []
	    for row in file:
	        d.append(row)
	    myfile.close()
	    firstRow=d[0]
	    return d




def detectSparsity(dataframe):
	i=0
	for row in dataframe:
		try:
			BTZ=[x if x=='' else float(x) for x in row[1:9]]
			LBH=[x if x=='' else float(x) for x in row[9:18]]
			LEN=[x if x=='' else float(x) for x in row[18:28]]
			NOT=[x if x=='' else float(x) for x in row[28:]]
			newBTZ=BTZ
			newLBH=LBH
			newLEN=LEN
			newNOT=NOT
			BTZnoGap=list(filter(lambda a: a != '', BTZ))
			LBHnoGap=list(filter(lambda b: b != '', LBH))
			LENnoGap=list(filter(lambda a: a != '', LEN))
			NOTnoGap=list(filter(lambda b: b != '', NOT))
			if all([BTZ.count('') >= 1, BTZ.count('') <5]):
				newBTZ=[x if x!='' else gmean(BTZnoGap) for x in BTZ]
			elif all([BTZ.count('') >=5, BTZ.count('') < 7]):
				#np.random.seed(200)
				try:
					k=2000
					centroids = {
						i+1: [np.random.randint(min(BTZnoGap),max(BTZnoGap)), np.random.randint(min(NOTnoGap),max(NOTnoGap))]
						for i in range(k)
					}
					newBTZ=[x if x!='' else list(centroids.values())[random.randint(0,50)][0] for x in BTZ]
				except ValueError:
					k=2000
					centroids = {
						i+1: [np.random.randint(2725302.821,12926007.57), np.random.randint(2423418.163,4006800.199)]
						for i in range(k)
					}

					newBTZ=[x if x!='' else list(centroids.values())[random.randint(0,50)][0] for x in BTZ]
			elif all([BTZ.count('') >= 7]):
				#random sampling from the min abundance of min of each 8 channels and max abundance of the min of each 8 channels
				#np.random.seed(200)
					k=2000
					centroids = {
						i+1: [np.random.randint(2725302.821,12926007.57), np.random.randint(2423418.163,4006800.199)]
						for i in range(k)
					}
					newBTZ=[x if x!='' else list(centroids.values())[random.randint(0,50)][0] for x in BTZ]
			dataframe[i][1:9]=newBTZ

			if all([LBH.count('') >= 1, LBH.count('') <5]):
				newLBH=[x if x!='' else gmean(LBHnoGap) for x in LBH]
			elif all([LBH.count('') >=5, LBH.count('') < 7]):
				#np.random.seed(200)
				try:
					k=2000
					centroids = {
						i+1: [np.random.randint(min(LBHnoGap),max(LBHnoGap)), np.random.randint(min(NOTnoGap),max(NOTnoGap))]
						for i in range(k)
					}

					newLBH=[x if x!='' else list(centroids.values())[random.randint(0,50)][0] for x in LBH]
				except ValueError:
					k=2000
					centroids = {
						i+1: [np.random.randint(2311917.816,3892664.915), np.random.randint(2423418.163,4006800.199)]
						for i in range(k)
					}

					newLBH=[x if x!='' else list(centroids.values())[random.randint(0,50)][0] for x in LBH]

			elif all([LBH.count('') >= 7]):
				#random sampling from the min abundance of min of each 8 channels and max abundance of the min of each 8 channels
				#np.random.seed(200)
				k=2000
				centroids = {
					i+1: [np.random.randint(2311917.816,3892664.915), np.random.randint(2423418.163,4006800.199)]
					for i in range(k)
				}
				newLBH=[x if x!='' else list(centroids.values())[random.randint(0,50)][0] for x in LBH]
			dataframe[i][9:18]=newLBH

			if all([LEN.count('') >= 1, LEN.count('') <5]):
				newLEN=[x if x!='' else gmean(LENnoGap) for x in LEN]
			elif all([LEN.count('') >=5, LEN.count('') < 8]):
				try:
					#np.random.seed(200)
					k=2000
					centroids = {
						i+1: [np.random.randint(min(LENnoGap),max(LENnoGap)), np.random.randint(min(NOTnoGap),max(NOTnoGap))]
						for i in range(k)
					}

					newLEN=[x if x!='' else list(centroids.values())[random.randint(0,50)][0] for x in LEN]
				except ValueError:
					k=2000
					centroids = {
						i+1: [np.random.randint(2724787.73,4638657.143), np.random.randint(2423418.163,4006800.199)]
						for i in range(k)
					}
					newLEN=[x if x!='' else list(centroids.values())[random.randint(0,50)][0] for x in LEN]
			elif all([LEN.count('') >= 8]):
				#random sampling from the min abundance of min of each 8 channels and max abundance of the min of each 8 channels
				#np.random.seed(200)
				k=2000
				centroids = {
					i+1: [np.random.randint(2724787.73,4638657.143), np.random.randint(2423418.163,4006800.199)]
					for i in range(k)
				}
				newLEN=[x if x!='' else list(centroids.values())[random.randint(0,50)][0] for x in LEN]
			dataframe[i][18:28]=newLEN

			if all([NOT.count('') >= 1, NOT.count('') <5]):
				newNOT=[x if x!='' else gmean(NOTnoGap) for x in NOT]
			elif all([NOT.count('') >=5, NOT.count('') < 8]):
				try:
					#np.random.seed(200)
					k=2000
					centroids = {
						i+1: [np.random.randint(min(NOTnoGap),max(NOTnoGap)), np.random.randint(min(NOTnoGap),max(NOTnoGap))]
						for i in range(k)
					}

					newNOT=[x if x!='' else list(centroids.values())[random.randint(0,50)][0] for x in NOT]
				except ValueError:
					#np.random.seed(200)
					k=2000
					centroids = {
						i+1: [np.random.randint(2423418.163,4006800.199), np.random.randint(2423418.163,4006800.199)]
						for i in range(k)
					}

					newNOT=[x if x!='' else list(centroids.values())[random.randint(0,50)][0] for x in NOT]
			elif all([NOT.count('') >= 8]):
				#random sampling from the min abundance of min of each 8 channels and max abundance of the min of each 8 channels
				#np.random.seed(200)
				k=2000
				centroids = {
					i+1: [np.random.randint(2423418.163,4006800.199), np.random.randint(2423418.163,4006800.199)]
					for i in range(k)
				}
				newNOT=[x if x!='' else list(centroids.values())[random.randint(0,50)][0] for x in NOT]
			dataframe[i][28:]=newNOT
	
		except TypeError:
			print('TypeError')
		i+=1

	return dataframe

def writeFile(infile, currData):
	fn = os.path.splitext(infile)[0] + 'FilledChannels_fixCluster.csv'
	with open(fn, 'w') as myfile:
	    outputFile = csv.writer(myfile)
	    outputFile.writerow(firstRow)
	    for site in currData:
	        outputFile.writerow(site)

def main():
	MetaData = readFile("CellMatrix_LEN_LBH_BTZ_NOT.csv")[1:]
	writeFile("CellMatrix_LEN_LBH_BTZ_NOT.csv", detectSparsity(MetaData))

if __name__=="__main__":
	main()
