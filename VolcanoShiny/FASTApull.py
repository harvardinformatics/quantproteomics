import sys
import csv
import os
import re
import numpy as np
from tkinter import *
import tkinter as tk

class ChannelFillApp(tk.Tk):
    def __init__(self):
        tk.Tk.__init__(self)
        self.v1=StringVar()
        self.v1.set("Enter FASTA filename")
        self.entry1 = tk.Entry(self, textvariable=self.v1)
        self.v2=StringVar()
        self.v2.set("Enter Accession filename")
        self.entry2 = tk.Entry(self, textvariable=self.v2)
        self.button = tk.Button(self, text="Assign", command=self.on_button)
        self.button.pack()
        self.button2 = tk.Button(self, text="Execute", command=self.destroy)
        self.button2.pack()
        self.entry1.pack()
        self.entry2.pack()
        
    def on_button(self):
    	global FASTAfile
    	global CSVfile

    	FASTAfile = self.entry1.get()
    	CSVfile = self.entry2.get()

    	print(self.entry1.get())
    	print(self.entry2.get())


app = ChannelFillApp()
app.mainloop()


FASTAfilename=FASTAfile
CSVfilename=CSVfile

def read_FASTA_strings(filename):

    with open(filename) as file:

        return file.read().split('>')[1:]

def read_FASTA_entries(filename):

    return [seq.partition('\n') for seq in read_FASTA_strings(filename)]

def read_FASTA_sequences(filename,faDict):
    for seq in read_FASTA_entries(filename):
    	try:
    		faDict[seq[0][0:].split('|')[1]]=seq[0]
    	except IndexError:
    		faDict[seq[0][0:]]=seq[0]
    return faDict
             


def readFile(infile):
	global firstRow
	fn = infile
	with open(fn, 'r') as myfile:
	    file = csv.reader(myfile, dialect='excel')
	    
	    d = []
	    for row in file:
	        d.append(row)
	    myfile.close()
	    return d

def matchGeneName(dat, faArr):
	newDat=[]
	for accession in dat:

		try:
			headerList=faArr.get(accession[0]).split(' ')
			matches = [x for x in headerList if 'GN=' in x]
			newDat.append([accession[0],matches[0].split('=')[1]])
		except:
			pass
	return newDat
	

def writeFile(infile, currData):
	fn = os.path.splitext(infile)[0] + 'OUTPUT.csv'
	with open(fn, 'w') as myfile:
	    outputFile = csv.writer(myfile)
	    for site in currData:
	        outputFile.writerow(site)

def main():
	humanDict={}
	fastaArray = read_FASTA_sequences(FASTAfilename,humanDict)
	da = readFile(CSVfilename)
	newAcc = matchGeneName(da, fastaArray)
	writeFile(CSVfilename, newAcc)
	print("your OUTPUT file has been written")

if __name__=="__main__":
	main()



