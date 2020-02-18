import sys
import csv
import os
import re
import numpy as np
import argparse
import statsmodels
from statsmodels.stats import multitest
from scipy.stats import ttest_ind
from scipy.stats import variation
from scipy.stats.mstats import gmean
from scipy.stats import rankdata
import pandas as pd
import random
import math
from tkinter import *
import tkinter as tk


class ChannelFillApp(tk.Tk):
    def __init__(self):
        tk.Tk.__init__(self)
        self.v1=StringVar()
        self.v1.set("Enter PSM filename")
        self.entry1 = tk.Entry(self, textvariable=self.v1)
        self.v2=StringVar()
        self.v2.set("Enter Replicate Number")
        self.entry2 = tk.Entry(self, textvariable=self.v2)
        self.v3=StringVar()
        self.v3.set("Enter Start Abundance Column")
        self.entry3 = tk.Entry(self, textvariable=self.v3)
        self.button = tk.Button(self, text="Assign", command=self.on_button)
        self.button.pack()
        self.button2 = tk.Button(self, text="Execute", command=self.destroy)
        self.button2.pack()
        self.entry1.pack()
        self.entry2.pack()
        self.entry3.pack()

    def on_button(self):
    	global PSMfile
    	global Replicate
    	global Column
    	PSMfile = self.entry1.get()
    	Replicate = int(self.entry2.get())*2
    	Column = self.entry3.get()

    	print(self.entry1.get())
    	print(self.entry2.get())
    	print(self.entry3.get())

app = ChannelFillApp()
app.mainloop()


begCol=int(Column)
midCol=int(int(Column)+(int(Replicate)/2))
endCol=int(Column)+int(Replicate)
reps=int(Replicate)


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

def removeBlankChannels(dfs):
	newdf=[]
	for row in dfs:
		#Change
		if row[begCol:endCol].count('')!=reps:
			newdf.append(row)
	return newdf


def detectSparsity(dataframe):
	i=0
	for row in dataframe:
		treatment=[x if x=='' else float(x) for x in row[begCol:midCol]]
		control=[x if x=='' else float(x) for x in row[midCol:endCol]]
		newTreatment=treatment
		newControl=control
		tnoGap=list(filter(lambda a: a != '', treatment))
		cnoGap=list(filter(lambda b: b != '', control))
		if all([treatment.count('') >= 1, treatment.count('') <reps, variation(tnoGap, axis = 0) > 0.3]):
			newTreatment=[x if x!='' else gmean(tnoGap) for x in treatment]
		elif all([treatment.count('') >= 1, treatment.count('') < reps-1, variation(tnoGap, axis = 0) <= 0.3]):
			#try:
			df = pd.DataFrame({'x':treatment, 'y':control})
			np.random.seed(200)
			k=5
			centroids = {
				i+1: [np.random.randint(min(tnoGap)-10,max(tnoGap)+10), np.random.randint(min(tnoGap)-10,max(tnoGap)+10)]
				for i in range(k)
			}
			print(centroids)
			newTreatment=[x if x!='' else list(centroids.values())[random.randint(0,4)][0] for x in treatment]
		elif all([treatment.count('') == reps, variation(tnoGap, axis = 0) == 0.0]):
			#try:
			df = pd.DataFrame({'x':treatment, 'y':control})
			np.random.seed(200)
			k=5
			centroids = {
				i+1: [np.random.randint(min(tnoGap)-25,max(tnoGap)+25), np.random.randint(min(tnoGap)-25,max(tnoGap)+25)]
				for i in range(k)
			}
			newTreatment=[x if x!='' else list(centroids.values())[random.randint(0,4)][0] for x in treatment]

		if all([control.count('') >= 1, control.count('') < reps, variation(cnoGap, axis = 0) > 0.3]):
			newControl=[x if x!='' else gmean(cnoGap) for x in control]
		elif all([control.count('') >= 1, control.count('') < reps-1, variation(cnoGap, axis = 0) <= 0.3]):
			#try:
			df = pd.DataFrame({'x':treatment, 'y':control})
			np.random.seed(200)
			k=5
			centroids = {
				i+1: [np.random.randint(min(cnoGap)-10,max(cnoGap)+10), np.random.randint(min(cnoGap)-10,max(cnoGap)+10)]
				for i in range(k)
			}
			print(centroids)
			newControl=[x if x!='' else list(centroids.values())[random.randint(0,4)][1] for x in control]
		elif all([control.count('') == 4, variation(cnoGap, axis = 0) == 0.0]):
			#try:
			df = pd.DataFrame({'x':treatment, 'y':control})
			np.random.seed(200)
			k=5
			centroids = {
				i+1: [np.random.randint(min(cnoGap)-25,max(cnoGap)+25), np.random.randint(min(cnoGap)-25,max(cnoGap)+25)]
				for i in range(k)
			}
			print(centroids)
			newControl=[x if x!='' else list(centroids.values())[random.randint(0,4)][1] for x in control]
		if all([treatment.count('') == reps, control.count('') == 0]):
			newControl=[100]*reps
			newTreatment=[1]*reps
		if all([treatment.count('') == 0, control.count('') == reps]):
			newControl=[1]*reps
			newTreatment=[100]*reps
		if treatment.count('') == reps:
			newTreatment=[1]*reps
		if control.count('') == reps:
			newControl=[1]*reps
		print('new values')
		print(newTreatment)
		print(newControl)
		if newControl!=[] and newTreatment!=[]:
			dataframe[i][begCol:midCol]=newTreatment
			dataframe[i][midCol:endCol]=newControl
			print('new values in data frame')
			print(dataframe[i][begCol:midCol])
			print(dataframe[i][midCol:endCol])
		i+=1
	return dataframe

def writeFile(infile, currData):
	fn = os.path.splitext(infile)[0] + 'FilledChannels.csv'
	with open(fn, 'w') as myfile:
	    outputFile = csv.writer(myfile)
	    outputFile.writerow(firstRow)
	    for site in currData:
	        outputFile.writerow(site)

def main():
	PSMdata = readFile(PSMfile)[1:]
	PSMfiltered = removeBlankChannels(PSMdata)
	writeFile(PSMfile, detectSparsity(PSMfiltered))

if __name__=="__main__":
	main()
