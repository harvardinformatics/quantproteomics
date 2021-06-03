#run on terminal with python3 IsotopeAnalysis.py 'C4 H7 N O4' 'negative' '[M-H]'

import sys
import csv
import os
import re
import numpy as np
import argparse

Infile=sys.argv[1]
mzFile=sys.argv[2]


def readFile(infile):
	fn = infile
	with open(fn, 'r') as myfile:
	    file = csv.reader(myfile, dialect='excel')
	    global firstRow
	    d = []
	    for row in file:
	        d.append(row)
	    myfile.close()
	    return d


def calcTopAbundance(mzdat, cmpdat):

	i=0
	highAbund=[]


	while i<len(mzdat[0]):
		
		mzLow=float(mzdat[0][i])-0.0005
		mzHigh=float(mzdat[0][i])+0.0005
		gatherIntensity=[]
		j=1
		while j<len(cmpdat):
			if mzLow<=float(cmpdat[j][1])<=mzHigh:
			
				gatherIntensity.append(float(cmpdat[j][2]))

				
			j+=1
		try:
			gatherIntensity.sort()
			highAbund.append(gatherIntensity[len(gatherIntensity)-1])


		except IndexError:
			highAbund.append(0)
		i+=1
	print(highAbund)
	return highAbund

def calcPerc(topAbund):
	sumAbund=sum(topAbund)
	percAr=[]
	for i in topAbund:
		percAr.append(round((100*(i/sumAbund)),1))
	return percAr



def writeFile(infile, topData, topPerc):
	fn = infile
	with open(fn, 'w', newline="") as myfile:
	    outputFile = csv.writer(myfile)
	    outputFile.writerow(['M','M+1','M+2','M+3','M+4'])
	    outputFile.writerow(topData)
	    outputFile.writerow(topPerc)

def main():

	compound_data=readFile(Infile)
	mzML_data=readFile(mzFile)
	topAbundances=calcTopAbundance(mzML_data,compound_data)
	percAbundances=calcPerc(topAbundances)
	writeFile('percentage_sum.csv',topAbundances,percAbundances)
	
if __name__=="__main__":
	main()

