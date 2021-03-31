#run on terminal with python3 IsotopeAnalysis.py 'C4 H7 N O4' 'negative' '[M-H]'

import sys
import csv
import os
import re
import numpy as np
import argparse

print(sys.argv)
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


def calcAbundancePerc(mzdat, cmpdat):
	i=0
	highAbund=[]
	while i<len(mzdat):
		#mzLow=float(mzdat[i])-5*(float(mzdat[i])/1000000)
		#mzHigh=float(mzdat[i])+5*(float(mzdat[i])/1000000)
		mzLow=132.0295698
		mzHigh=132.0308902
		#print('mzLow:'+str(mzLow))
		#print('mzHigh:'+str(mzHigh))
		gatherIntensity=[]
		j=1
		while j<len(cmpdat):
			#print(float(cmpdat[j][0]))
			if mzLow<=float(cmpdat[j][0])<=mzHigh:
				gatherIntensity.append(float(cmpdat[j][1]))
				#print(float(cmpdat[j][1]))
			
			j+=1
		gatherIntensity.sort()
		highAbund.append(gatherIntensity[len(gatherIntensity)-1])

		i+=1
	return highAbund




def writeFile(infile, currData):
	fn = os.path.splitext(infile)[0] + '_.csv'
	with open(fn, 'w', newline="") as myfile:
	    outputFile = csv.writer(myfile)
	    outputFile.writerow(firstRow)
	    for site in currData:
	        outputFile.writerow(site)

def main():

	compound_data=readFile(Infile)
	mzML_data=readFile(mzFile)[0]
	print(calcAbundancePerc(mzML_data,compound_data))
	
if __name__=="__main__":
	main()

