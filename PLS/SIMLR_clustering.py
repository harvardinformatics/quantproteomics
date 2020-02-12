import csv
from pandas import DataFrame
import statsmodels.api as sm 
from statsmodels.api import graphics
from sklearn import datasets
import matplotlib.pyplot as plt 
import math
import seaborn as sns
from statsmodels.formula.api import ols


def readFile(infile):
	global firstRow
	fn = infile
	with open(fn, 'r') as myfile:
	    file = csv.reader(myfile, dialect='excel')
	    x = []
	    y = []
	    z = []
	    for row in file:
	        x.append(math.log2(float(row[7])))
	        y.append(math.log2(float(row[18])))
	        z.append(str(row[2]))

	    myfile.close()
	    return x, y, z

X, y, z  = readFile('191213_Myel_LEN_C1_20minBOTH_LEN_Matched.csv')

metabolomics={'NOT':X, 'LBH':y, 'RT':z}
df=DataFrame(metabolomics, columns = ['NOT','LBH','RT'])
print(df)


model = sm.OLS(y, X).fit()
predictions = model.predict(X) # make the predictions by the model

# Print out the statistics
print(model.summary())
fig, ax = plt.subplots(figsize=(8,6))
fig = sm.graphics.influence_plot(model, ax=ax)
plt.show()

fig, ax = plt.subplots(figsize=(8,6))
fig = sm.graphics.plot_leverage_resid2(model, ax=ax)
plt.show()

fig, ax = plt.subplots(figsize=(8,6))
fig = sm.graphics.plot_partregress("NOT", "LBH", ["RT"], data=df, ax=ax)
plt.show()




