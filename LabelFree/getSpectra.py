from pyopenms import *
import csv



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
	x=1
	while x < len(inDat):
		PSM=inDat[x]
		exp=MSExperiment()
		fn=PSM[1].split('.')[0]+'.mzML'
		MzMLFile().load(fn, exp)
		gAr=[]
		sumAr=[]
		for spectrum in exp:
		    for peak in spectrum:
		    	if all([spectrum.getRT() >=(float(PSM[2])-0.005)*60, spectrum.getRT() <=(float(PSM[2])+0.005)*60]):
		        	gAr.append([spectrum.getRT(), peak.getMZ(), peak.getIntensity(),spectrum.getRT()/60])

		        	sumAr.append(peak.getIntensity())
		inDat[x].append(sum(sumAr))
		print(sum(sumAr))
		x+=1
	return inDat

def writeFile(dframe):
	fn = 'test.csv'
	with open(fn, 'w', newline="") as myfile:
	    outputFile = csv.writer(myfile)
	    for i in dframe:
	    	outputFile.writerows([i])



def main():
	
	psm_data=readFile('JunB_Testing.csv')
	sumAbundances=getAreas(psm_data)
	
	writeFile(sumAbundances)
	
if __name__=="__main__":
	main()