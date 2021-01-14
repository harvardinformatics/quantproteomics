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
from itertools import permutations




def readFile(infile):
	fn = infile
	with open(fn, 'r') as myfile:
	    file = csv.reader(myfile, dialect='excel')
	    d = []
	    for row in file:
	    	if 'Confidence' in row:
	    		pass
	    	else:
	    		d.append(row[29:37])
	    myfile.close()
	return d


def createMatrix(dataM):
	col=len(dataM)-1
	row=len(dataM[0])-1
	perc2Data=int(col*row*0.02)
	perc5Data=int(col*row*0.05)
	matrixStorage={-1:dataM}
	dataMnew=dataM
	numBlanks=random.randint(perc2Data,perc5Data)
	
	for j in range(numBlanks):
		cRand=random.randint(0, col)
		rRand=random.randint(0, row)
		dataMnew[cRand][rRand]=''
	

	return dataMnew

def createHashtable(dataNewm):
	hasht={}
	for i in range(100):
		hasht[i]=dataNewm
	return hasht

def blankPredictor(obsMatrix, hashMatrix):
	percMatrix=obsMatrix
	i=0
	while i < len(obsMatrix):
		j=0
		blank=0
		nonblank=0
		while j < len(obsMatrix[i]):
			for key,value in hashMatrix.items():
				
				
				if value[i][j] == '':
					blank+=1
				else:
					nonblank+=1
				percMatrix[i][j]=100*blank/(blank+nonblank)
			j+=1
		i+=1
	return percMatrix


def writeFile(infile, currData):
	fn = os.path.splitext(infile)[0] + '_Imputed.csv'
	with open(fn, 'w', newline="") as myfile:
	    outputFile = csv.writer(myfile)
	    outputFile.writerow(firstRow)
	    for site in currData:
	        outputFile.writerow(site)

def main():
	PSMmatrix = createHashtable(createMatrix(PSMdata))



if __name__=="__main__":
	main()
