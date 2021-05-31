#run on terminal with python3 IsotopeAnalysis.py 'C4 H7 N O4' 'negative' '[M-H]'

import sys
import csv
import os
import re
import numpy as np
import argparse


Inpsm=sys.argv[1]


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


def calcTopAbundance(mzdat):
	j=1
	while j<len(mzdat):

		gatherIntensity=[]
		mzLow=float(mzdat[j][16])-10
		mzHigh=float(mzdat[j][16])+10
		
		for Infile in os.listdir('./spectra'):
			try:
				mzmlFile=Infile.split(' ')[2]+'.raw'
				scanFile=Infile.split(' ')[3]
			except IndexError:
				pass

			if all([mzdat[j][26] == scanFile, mzdat[j][27] == mzmlFile]):
				try:
					compound_data=readFile('./spectra/'+Infile)
					k=1
					while k<len(compound_data):
						if mzLow<=float(compound_data[k][0])<=mzHigh:
							gatherIntensity.append(float(compound_data[j][1]))
						k+=1
				except UnicodeDecodeError:
					pass

				
		try:
			gatherIntensity.sort()
			mzdat[j].append(sum(gatherIntensity))
		except IndexError:
			mzdat[j].append('NA')
		j+=1


	return mzdat




def writeFile(infile, topData):
	
	fn = infile
	with open(fn, 'w', newline="") as myfile:
	    outputFile = csv.writer(myfile)
	    outputFile.writerows(topData)
	

def main():
	
	psm_data=readFile(Inpsm)
	topAbundances=calcTopAbundance(psm_data)
	
	writeFile('peakAreas.csv',topAbundances)
	
if __name__=="__main__":
	main()

