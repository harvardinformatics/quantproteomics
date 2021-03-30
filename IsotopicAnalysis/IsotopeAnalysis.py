#run on terminal with python3 IsotopeAnalysis.py 'C4 H7 N O4'

import sys

Compound=sys.argv[1]

def createEleDict():
	elementDict={'H':1.00797,'He':4.00260,'Li':6.941,'Be':9.01218,'B':10.81,'C':12.011,'N':14.0067,'O':15.9994,'F':18.998403,'Ne':20.179,
	'Ne':20.179,'Na':22.98977,'Mg':24.305,'Al':26.98154,'Si':28.0855,'P':30.97376,'S':32.06,'Cl':35.453,'K':39.0983,'Ar':39.948,'Ca':40.08,
	'Sc':44.9559,'Ti':47.90,'V':50.9415,'Cr':51.996,'Mn':54.9380,'Fe':55.847,'Ni':58.70,'Co':58.9332,'Cu':63.546,'Zn':65.38,'Ga':69.72,
	'Ge':72.59,'As':74.9,'Se':78.96,'Br':79.904,'Kr':83.80,'Rb':85.4678,'Sr':87.62,'Y':88.9059,'Zr':91.22,'Nb':92.9064,'Mo':95.94,
	'Tc':98,'Ru':101.07,'Rh':102.9055,'Pd':106.4,'Ag':107.868,'Cd':112.41,'In':114.82,'Sn':118.69,'Sb':121.75,'I':126.9045,'Te':127.60}
	return elementDict


def createIsotopes(cmp, eleDict):
	xC=cmp.split(' ')
	isotopes=[]
	C=0
	mz=0
	for element in xC:
		
		for x in element:
			
			if x.isnumeric():
				ele=element[:element.find(str(x))]
				eleNum=int(element[element.find(str(x)):])
				mz+=eleDict[ele]*eleNum

				
				if 'C' in element:
					C=eleNum
	isotopes.append(mz)
	i=0
	Mplus=mz
	while i<C:
		Mplus+=1.00336
		isotopes.append(Mplus)                       
		i+=1
	return isotopes

def main():
	print(createIsotopes(Compound, createEleDict()))

if __name__=="__main__":
	main()