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
from dupl_matrix import createMatrix, createHashtable, blankPredictor, readFile

from missingpy import MissForest


begCol=int(sys.argv[3])
midCol=int(int(sys.argv[3])+(int(sys.argv[2])))
endCol=midCol+int(sys.argv[4])
reps=int(sys.argv[2])
reps2=int(sys.argv[4])



def readFileNAN(infile):
	print('Your file is being processed, stay tuned!')
	global firstRow
	fn = infile
	with open(fn, 'r') as myfile:
	    file = csv.reader(myfile, dialect='excel')
	    d = []
	    annotation = []
	    for row in file:
	    	if row[begCol:endCol].count('')==reps+reps2:
	    		pass
	    	else:
		    	if 'Confidence' in row:
		    		firstRow = row
		    	else:
		    		i=begCol
		    		while i < endCol+1:

		    			try:
		    				row[i]=float(row[i])
		    			except:
		    				row[i]=np.nan
		    			
		    			i+=1
		    		d.append(row[begCol:endCol])
		    		annotation.append(row)
	    myfile.close()
	return d, annotation


def imputeMatrix(dataM):
	
	nan=np.nan
	imputer = MissForest()
	dataT = imputer.fit_transform(dataM)
	return dataT

def blankDecider(pmatrix, bmatrix, imatrix):
	x=0

	nan=np.nan
	while x<len(pmatrix):
		y=0
		
		while y<len(pmatrix[x]):
			if all([pmatrix[x][y] is np.nan, bmatrix[x][y]>=0.6]):

		
				imatrix[x,y]=np.nan

			y+=1
		x+=1
	return imatrix


def writeFile(infile, currData, annotData):

	cData = currData.tolist()
	fn = os.path.splitext(infile)[0] + '_Imputed_missForest.csv'
	with open(fn, 'w', newline="") as myfile:
	    outputFile = csv.writer(myfile)
	    outputFile.writerow(firstRow)
	    x=0
	    while x<len(cData):
	        outputFile.writerow(annotData[x][0:begCol]+cData[x]+annotData[x][endCol+1:])
	        x+=1
	print('Your file has been output with appropriate imputation of missing values.')




def main():
	PSMdataNAN, PSMannot = readFileNAN('200918A_DMSAM07342_DKMIVB20-PSM.csv')
	PSMdata, PSMannot = readFileNAN('200918A_DMSAM07342_DKMIVB20-PSM.csv')
	impMatrix = imputeMatrix(PSMdata)
	PSMmatrix = createHashtable(createMatrix(PSMdata))
	blankMatrix = blankPredictor(PSMdata,PSMmatrix)
	
	writeFile('200918A_DMSAM07342_DKMIVB20-PSM.csv',blankDecider(PSMdataNAN,blankMatrix,impMatrix), PSMannot)



if __name__=="__main__":
	main()
