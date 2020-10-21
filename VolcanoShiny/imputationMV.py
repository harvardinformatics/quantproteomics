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
from collections import defaultdict




begCol=int(sys.argv[3])
midCol=int(int(sys.argv[3])+(int(sys.argv[2])))
endCol=midCol+int(sys.argv[5])
reps=int(sys.argv[2])
reps2=int(sys.argv[5])

level=int(sys.argv[4])


def readFile(infile):
	print('Your file is being processed, stay tuned!')
	global firstRow
	global noiseTmin
	global noiseTmax
	global noiseCmin
	global noiseCmax
	fn = infile
	columns = defaultdict(list)
	with open(fn, 'r') as myfile:
	    file = csv.reader(myfile, dialect='excel')
	    
	    d = []
	    for row in file:
	        d.append(row)
	        for (i,v) in enumerate(row):
	        	columns[i].append(v)
	    myfile.close()
	    
	    firstRow=d[0]
	    x=begCol
	    noiseTreatment=[]

	    while x<midCol:
	    	columns[x]=list(filter(lambda a: a != '', columns[x][1:]))
	    	columns[x]=[float(i) for i in columns[x]]
	    	#print(columns[x])
	    	noiseTreatment.append(min(columns[x]))
	    	x+=1
	    y=midCol
	    noiseControl=[]
	    while y<endCol:
	    	columns[y]=list(filter(lambda a: a != '', columns[y][1:]))
	    	columns[y]=[float(i) for i in columns[y]]
	    	noiseControl.append(min(columns[y]))
	    	y+=1
	    noiseTmin=min(noiseTreatment)
	    noiseTmax=max(noiseTreatment)
	    noiseCmin=min(noiseControl)
	    noiseCmax=max(noiseControl)

	    return d

def removeBlankChannels(dfs):
	newdf=[]
	for row in dfs:
		#Change
		if row[begCol:endCol].count('')!=reps+reps2:
			newdf.append(row)
	return newdf


def detectSparsity(dataframe):
	
	i=0
	for row in dataframe:
		treatment=[x if x=='' else float(x) for x in row[begCol:midCol]]
		control=[x if x=='' else float(x) for x in row[midCol:endCol]]
		newTreatment=treatment
		newControl=control
		tnoGap=list(filter(lambda a: a != '', treatment))
		cnoGap=list(filter(lambda b: b != '', control))
		if treatment.count('') <reps-2:
			newTreatment=[x if x!='' else gmean(tnoGap) for x in treatment]
		elif treatment.count('') < reps-1:
			try:
				k=2000
				centroids = {
					i+1: [np.random.randint(min(tnoGap),max(tnoGap)), np.random.randint(min(cnoGap),max(cnoGap))]
					for i in range(k)
				}

				newTreatment=[x if x!='' else list(centroids.values())[random.randint(0,50)][0] for x in treatment]

			except ValueError:
				if level == 0:
					pass
				elif level == 1:
					if variation(tnoGap, axis = 0) > 0.3:
						k=2000
						centroids = {
							i+1: [np.random.randint(noiseTmin,noiseTmax), np.random.randint(noiseCmin, noiseCmax)]
							for i in range(k)
						}
						#print([MORnoGap,CTLnoGap])
						#print('centroids')
						#print(centroids)
						newTreatment=[x if x!='' else list(centroids.values())[random.randint(0,50)][0] for x in treatment]


				else:
					k=2000
					centroids = {
						i+1: [np.random.randint(noiseTmin,noiseTmax), np.random.randint(noiseCmin, noiseCmax)]
						for i in range(k)
					}
					#print([MORnoGap,CTLnoGap])
					#print('centroids')
					#print(centroids)
					newTreatment=[x if x!='' else list(centroids.values())[random.randint(0,50)][0] for x in treatment]



		elif treatment.count('') <= reps:
			if any([level == 0, level == 1]):
				pass
			else:
				
				k=2000
				centroids = {
					i+1: [np.random.randint(noiseTmin,noiseTmax), np.random.randint(noiseCmin, noiseCmax)]
					for i in range(k)
				}
				#print([MORnoGap,CTLnoGap])
				#print('centroids')
				#print(centroids)
				newTreatment=[x if x!='' else list(centroids.values())[random.randint(0,50)][0] for x in treatment]



		if control.count('') < reps2-2:
			newControl=[x if x!='' else gmean(cnoGap) for x in control]
		elif control.count('') < reps2-1:
			try:
				k=2000
				centroids = {
					i+1: [np.random.randint(min(tnoGap),max(tnoGap)), np.random.randint(min(cnoGap),max(cnoGap))]
					for i in range(k)
				}

				newControl=[x if x!='' else list(centroids.values())[random.randint(0,50)][1] for x in control]

			except ValueError:
				if level == 0:
					pass
				elif level == 1:
					if variation(cnoGap, axis = 0) > 0.3:
						k=2000
						centroids = {
							i+1: [np.random.randint(noiseTmin,noiseTmax), np.random.randint(noiseCmin, noiseCmax)]
							for i in range(k)
						}
						#print([MORnoGap,CTLnoGap])
						#print('centroids')
						#print(centroids)
						newControl=[x if x!='' else list(centroids.values())[random.randint(0,50)][1] for x in control]

				else:
					k=2000
					centroids = {
						i+1: [np.random.randint(noiseTmin,noiseTmax), np.random.randint(noiseCmin, noiseCmax)]
						for i in range(k)
					}
					#print([MORnoGap,CTLnoGap])
					#print('centroids')
					#print(centroids)
					newControl=[x if x!='' else list(centroids.values())[random.randint(0,50)][1] for x in control]

		elif control.count('') <= reps2:
			if any([level == 0, level == 1]):
				pass
			else:
				
				k=2000
				centroids = {
					i+1: [np.random.randint(noiseTmin,noiseTmax), np.random.randint(noiseCmin, noiseCmax)]
					for i in range(k)
				}
				#print([MORnoGap,CTLnoGap])
				#print('centroids')
				#print(centroids)
				newControl=[x if x!='' else list(centroids.values())[random.randint(0,50)][1] for x in control]


		if newControl!=[] and newTreatment!=[]:
			dataframe[i][begCol:midCol]=newTreatment
			dataframe[i][midCol:endCol]=newControl

		i+=1
	return dataframe

def writeFile(infile, currData):
	fn = os.path.splitext(infile)[0] + '_Imputed.csv'
	with open(fn, 'w', newline="") as myfile:
	    outputFile = csv.writer(myfile)
	    outputFile.writerow(firstRow)
	    for site in currData:
	        outputFile.writerow(site)
	print('Your file has been output with appropriate imputation of missing values.')

def main():

	PSMdata = readFile(sys.argv[1])[1:]
	PSMfiltered = removeBlankChannels(PSMdata)
	writeFile(sys.argv[1], detectSparsity(PSMfiltered))


if __name__=="__main__":
	main()
