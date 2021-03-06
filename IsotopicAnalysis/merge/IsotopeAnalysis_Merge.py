#run on terminal with python3 IsotopeAnalysis.py 'C4 H7 N O4' 'negative' '[M-H]'
from pyopenms import *
import sys
import csv

Compound=sys.argv[1]
IonMode=sys.argv[2]
AdductAddition=sys.argv[3]

mzFile=sys.argv[4]
rTime=float(sys.argv[5])


Infile='compound_spectra.csv'
mzFileCAP='mz_calc.csv'

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


def createEleDict():
	elementDict={'H':1.007825,'He':4.00260,'Li':6.941,'Be':9.01218,'B':10.81,'C':12,'N':14.003074,'O':15.9949146,'F':18.998403,'Ne':20.179,
	'Ne':20.179,'Na':22.98977,'Mg':24.305,'Al':26.98154,'Si':28.0855,'P':30.97376,'S':32.06,'Cl':35.453,'K':39.0983,'Ar':39.948,'Ca':40.08,
	'Sc':44.9559,'Ti':47.90,'V':50.9415,'Cr':51.996,'Mn':54.9380,'Fe':55.847,'Ni':58.70,'Co':58.9332,'Cu':63.546,'Zn':65.38,'Ga':69.72,
	'Ge':72.59,'As':74.9,'Se':78.96,'Br':79.904,'Kr':83.80,'Rb':85.4678,'Sr':87.62,'Y':88.9059,'Zr':91.22,'Nb':92.9064,'Mo':95.94,
	'Tc':98,'Ru':101.07,'Rh':102.9055,'Pd':106.4,'Ag':107.868,'Cd':112.41,'In':114.82,'Sn':118.69,'Sb':121.75,'I':126.9045,'Te':127.60}
	# otherDict={'Xe':131.30,'Cs':132.9054,'Ba':137.33,'La':138.9055,'Ce':140:12,'Pr':140.9077,'Nd':144.24,'Sm':150.4,'Eu':151.96,'Gd':151.96,
	# 'Td':158.9254,'Dy':162.50,'Ho':164.9304,'Er':167.26,'Tm':168.9342,'Yb':173.04,'Lu':174.967,'Hf':178.49,'Ta':180.9479,'W':183.85,'Re':186.207,
	# 'Os':190.2,'Ir':192.22,'Pt':195.09,'Au':196.9665,'Hg':200.59,'Tl':204.37,'Pb':207.2,'Bi':208.9804,'Po':209.0,'At':210.0,'Rn':222.0,'Fr':223.0,
	# 'Ra':226.0254,'Ac':227.0278,'Pa':231.0359,'Th':232.0381,'Np':237.0482,'U':238.029}
	return elementDict

def createPosAddict():
	def calcMplusH(M):
		return M+1.007825-0.0005485799
	def calcMplusNa(M):
		return M+22.989218
	posAddict={'[M+H]':calcMplusH, '[M+Na]':calcMplusNa}
	return posAddict

def createNegAddict():
	def calcMminusH(M):
		return M-1.007825+0.0005485799

	negAddict={'[M-H]':calcMminusH}
	return negAddict


def createIsotopes(cmp, eleDict):
	xC=cmp.split('_')
	isotopes=[]
	C=0
	mz=0
	for element in xC:
		
		for x in element:
			
			if x.isnumeric():
				ele=element[:element.find(str(x))]
				eleNum=int(element[element.find(str(x)):])
				print(eleDict[ele])
				print(eleNum)
				mz+=eleDict[ele]*eleNum
				print('mz: '+str(mz))
				
				if 'C' in element:
					C=eleNum
			elif element.isalpha():
				ele=x
				eleNum=1
				print(eleDict[ele])
				print(eleNum)
				mz+=eleDict[ele]*eleNum
				print('mz: '+str(mz))

	isotopes.append(mz)
	i=0
	Mplus=mz
	while i<C:
		Mplus+=1.0033548
		isotopes.append(Mplus)                       
		i+=1
	return isotopes

def adductCalc(isoAr, ionMode, addct):
	newIsoDist=[]
	for iso in isoAr:

		if ionMode=='negative':

			newIsoDist.append(createNegAddict()[addct](iso))
		elif ionMode=='positive':

			newIsoDist.append(createPosAddict()[addct](iso))
	return newIsoDist


def getAreas(inDat):
	exp=MSExperiment()
	MzMLFile().load(inDat, exp)
	gAr=[]
	for spectrum in exp:
	    for peak in spectrum:
	    	if all([spectrum.getRT() >=(rTime-0.025)*60, spectrum.getRT() <=(rTime+0.025)*60]):
	        	gAr.append([spectrum.getRT(), peak.getMZ(), peak.getIntensity(),spectrum.getRT()/60])
	return gAr

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


def writeFileGS(infile,dframe):
	fn = infile
	with open(fn, 'w', newline="") as myfile:
	    outputFile = csv.writer(myfile)
	    for i in dframe:
	    	outputFile.writerows([i])

def writeFileCAP(infile, topData, topPerc):
	fn = infile
	with open(fn, 'w', newline="") as myfile:
	    outputFile = csv.writer(myfile)
	    outputFile.writerow(['M','M+1','M+2','M+3','M+4'])
	    outputFile.writerow(topData)
	    outputFile.writerow(topPerc)


def writeFile(infile, currData):
	fn = infile
	with open(fn, 'w', newline="") as myfile:
	    outputFile = csv.writer(myfile)
	    outputFile.writerow(currData)

def main():
	isos=createIsotopes(Compound, createEleDict())
	writeFile('mz_calc.csv',adductCalc(isos,IonMode,AdductAddition))
	sumAbundances=getAreas(mzFile)
	writeFileGS('compound_spectra.csv', sumAbundances)
	compound_data=readFile(Infile)
	mzML_data=readFile(mzFileCAP)
	topAbundances=calcTopAbundance(mzML_data,compound_data)
	percAbundances=calcPerc(topAbundances)
	writeFileCAP('percentage_sum.csv',topAbundances,percAbundances)

if __name__=="__main__":
	main()