from pyopenms import *
import csv
import sys

mzFile=sys.argv[1]
rTime=float(sys.argv[2])


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


def getAreas(inDat):
	exp=MSExperiment()
	MzMLFile().load(inDat, exp)
	gAr=[]
	for spectrum in exp:
	    for peak in spectrum:
	    	if all([spectrum.getRT() >=(rTime-0.025)*60, spectrum.getRT() <=(rTime+0.025)*60]):
	        	gAr.append([spectrum.getRT(), peak.getMZ(), peak.getIntensity(),spectrum.getRT()/60])
	return gAr

def writeFile(infile,dframe):
	fn = infile
	with open(fn, 'w', newline="") as myfile:
	    outputFile = csv.writer(myfile)
	    for i in dframe:
	    	outputFile.writerows([i])


def main():
	sumAbundances=getAreas(mzFile)
	writeFile('compound_spectra.csv', sumAbundances)

if __name__=="__main__":
	main()