import sys
import csv
import os
import re
import numpy as np
import argparse

psmfile=sys.argv[1]
protfile=sys.argv[2]
psmAccIndex=int(sys.argv[3])
protAccIndex=int(sys.argv[4])

def readFilePSM(infile):
	fn = infile
	with open(fn, 'r') as myfile:
	    file = csv.reader(myfile, dialect='excel')
	    global firstRow
	    d = []
	    for row in file:
	        d.append(row)
	    myfile.close()
	    firstRow = d[0]
	    return d

def readFileProtein(infile):
	fn = infile
	with open(fn, 'r') as myfile:
	    file = csv.reader(myfile, dialect='excel')
	    
	    d = []
	    for row in file:
	        d.append(row)
	    myfile.close()
	    return d

def match(dat, pd, psmAcc, protAcc):
	filterdat=[]
	i=0
	while i<len(dat):
		j=0
		while j<len(pd):
			if dat[i][psmAcc] == pd[j][protAcc]:
				filterdat.append(dat[i])
				break


			j+=1
		i+=1
	return filterdat




def writeFile(infile, currData):
	fn = os.path.splitext(infile)[0] + '_PDfilters.csv'
	with open(fn, 'w', newline="") as myfile:
	    outputFile = csv.writer(myfile)
	    #outputFile.writerow(firstRow)
	    for site in currData:
	        outputFile.writerow(site)

def main():
	protDesc=readFileProtein(protfile)
	data=readFilePSM(psmfile)
	
	writeFile(psmfile, match(data, protDesc, psmAccIndex, protAccIndex))
	print('Your PSM file has been output with appropriate PD filters.')

if __name__=="__main__":
	main()

