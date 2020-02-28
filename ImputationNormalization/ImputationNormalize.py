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


def removeBlankChannels(dfs):
	newdf=[]
	for row in dfs:
		#Change
		if row[15:23].count('')!=8:
			newdf.append(row)
	return newdf


def detectSparsity(dataframe):
	i=0
	for row in dataframe:
		try:
			MOR=[x if x=='' else float(x) for x in row[15:19]]
			CTL=[x if x=='' else float(x) for x in row[19:23]]

			newMOR=MOR
			newCTL=CTL
		
			MORnoGap=list(filter(lambda a: a != '', MOR))
			CTLnoGap=list(filter(lambda b: b != '', CTL))

			if MOR.count('') == 1:
				newMOR=[x if x!='' else gmean(MORnoGap) for x in MOR]
			elif MOR.count('') ==2:
				#np.random.seed(200)
				try:
					k=2000
					centroids = {
						i+1: [np.random.randint(min(MORnoGap),max(MORnoGap)), np.random.randint(min(CTLnoGap),max(CTLnoGap))]
						for i in range(k)
					}
					#print([MORnoGap,CTLnoGap])
					#print('centroids')
					#print(centroids)
					newMOR=[x if x!='' else list(centroids.values())[random.randint(0,50)][0] for x in MOR]

				except ValueError:
					k=2000
					centroids = {
						i+1: [np.random.randint(485.5877075,529.3711548), np.random.randint(472.7730103,558.460144)]
						for i in range(k)
					}
					#print([MORnoGap,CTLnoGap])
					#print('centroids')
					#print(centroids)
					newMOR=[x if x!='' else list(centroids.values())[random.randint(0,50)][0] for x in MOR]
			elif MOR.count('') >= 3:
				pass
				#random sampling from the min abundance of min of each 8 channels and max abundance of the min of each 8 channels
				#np.random.seed(200)
					# k=2000
					# centroids = {
					# 	i+1: [np.random.randint(485.5877075,529.3711548), np.random.randint(472.7730103,558.460144)]
					# 	for i in range(k)
					# }
					# #print([MORnoGap,CTLnoGap])
					# #print('centroids')
					# #print(centroids)
					# newMOR=[x if x!='' else list(centroids.values())[random.randint(0,50)][0] for x in MOR]
			dataframe[i][15:19]=newMOR

			if CTL.count('') == 1:
				newCTL=[x if x!='' else gmean(CTLnoGap) for x in CTL]
			elif CTL.count('') ==2:
				#np.random.seed(200)
				try:
					k=2000
					centroids = {
						i+1: [np.random.randint(min(MORnoGap),max(MORnoGap)), np.random.randint(min(CTLnoGap),max(CTLnoGap))]
						for i in range(k)
					}
					#print([MORnoGap,CTLnoGap])
					#print('centroids')
					#print(centroids)
					newCTL=[x if x!='' else list(centroids.values())[random.randint(0,50)][1] for x in CTL]
				except ValueError:
					k=2000
					centroids = {
						i+1: [np.random.randint(485.5877075,529.3711548), np.random.randint(472.7730103,558.460144)]
						for i in range(k)
					}
					#print([MORnoGap,CTLnoGap])
					#print('centroids')
					#print(centroids)
					newCTL=[x if x!='' else list(centroids.values())[random.randint(0,50)][1] for x in CTL]
			elif CTL.count('') >= 3:
				pass
				#random sampling from the min abundance of min of each 8 channels and max abundance of the min of each 8 channels
				#np.random.seed(200)
					# k=2000
					# centroids = {
					# 	i+1: [np.random.randint(485.5877075,529.3711548), np.random.randint(472.7730103,558.460144)]
					# 	for i in range(k)
					# }
					# #print([MORnoGap,CTLnoGap])
					# #print('centroids')
					# #print(centroids)
					# newCTL=[x if x!='' else list(centroids.values())[random.randint(0,50)][1] for x in CTL]
			dataframe[i][19:23]=newCTL

	
		except TypeError:
			print('TypeError')
		i+=1

	return dataframe

def writeFile(infile, currData):
	fn = os.path.splitext(infile)[0] + '_ImputedHalf.csv'
	with open(fn, 'w') as myfile:
	    outputFile = csv.writer(myfile)
	    outputFile.writerow(firstRow)
	    for site in currData:
	        outputFile.writerow(site)

def main():
	MetaData = readFile("200113L_SAM6526_BY53_CC_PO4_pH_All_KB_Percolator_sitelevel.csv")[1:]
	xData=removeBlankChannels(MetaData)
	writeFile("200113L_SAM6526_BY53_CC_PO4_pH_All_KB_Percolator_sitelevel.csv", detectSparsity(xData))

if __name__=="__main__":
	main()
