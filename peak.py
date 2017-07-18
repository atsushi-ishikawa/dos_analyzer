import numpy as np
import peakutils 
from peakutils.plot import plot as pplot
import matplotlib.pylab as plt

x = []
y = []

rdf_file = "asdf.txt"

# data = np.genfromtxt(rdf_file, delimiter="\s")
data = np.loadtxt(rdf_file)

x = data[:,0]
y = data[:,1]

indexes = peakutils.indexes(y, thres=0.01, min_dist=1)
pplot(x,y,indexes)
plt.show()

print(x[indexes],y[indexes])
