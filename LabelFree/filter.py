import sys
import csv
import os
import re
import numpy as np
import argparse



def readFile(infile):
	global firstRow
	fn = infile
	with open(fn, 'r') as myfile:
	    file = csv.reader(myfile, dialect='excel')
	    
	    d = []
	    for row in file:
	        d.append(row)
	    myfile.close()
	    firstRow=d[0]
	    return d


def countPSMs(PSMlist, PSMseq):
	countsPSM=0
	i=0
	while i<len(PSMlist):
		if PSMseq == PSMlist[i]:
			countsPSM+=1
		i+=1
	return countsPSM


def PSMposition(modi,prot):
	modif=modi.split(';')
	if len(modif) == 1:
		modPos=[]
		modifications=modif[0].split('(')[0][1:]

		modPos.append(int(modifications)+int(prot)-1)
	else:
		modPos=[]
		for mod in modif:
			modifications=mod.strip(' ').split('(')[0][1:]

			modPos.append(int(modifications)+int(prot)-1)
	return modPos, len(modif)


def filterPSMs(dat):
	i=1
	while i<len(dat):
		if any(['Deamidated' in dat[i][5], 'Oxidation' in dat[i][5]]):
			pass

		if countPSMs(dat[4],dat[i][4]) <= 2:
			pass

		dat[i].append(PSMposition(dat[i][5],dat[i][44])[0])
		dat[i].append(PSMposition(dat[i][5],dat[i][44])[1])
		i+=1
	return dat




def writeFile(infile, currData):
	fn = os.path.splitext(infile)[0] + '_OUTPUT.csv'
	fr=firstRow+['PSM position','Glycosite Number']
	cd=currData.pop(0)
	with open(fn, 'w') as myfile:
	    outputFile = csv.writer(myfile)
	    outputFile.writerow(fr)
	    for site in currData:
	    	outputFile.writerow(site)

def main():
	data=readFile('201123L_SAM07668_BY90_DR294_DR294A_CF_final_JUNB_PSMs.csv')

	writeFile('201123L_SAM07668_BY90_DR294_DR294A_CF_final_JUNB_PSMs.csv', filterPSMs(data))

if __name__=="__main__":
	main()

