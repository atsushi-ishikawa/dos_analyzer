from ase.calculators.vasp import VaspDos
import sys
import matplotlib.pylab as plt
import numpy as np

###
def gaussian(x,x0,a,b):
	"""
		y = a*exp(-b*(x-x0)**2)
	"""
	import numpy as np
	x = np.array(x)
	y = np.exp( -b*(x-x0)**2 )
	y = a*y
	return y
###

argvs  = sys.argv
doscar = argvs[1]

dos = VaspDos(doscar=doscar)

x = dos.energy
y = dos.dos

len = len(x)
y2  = np.zeros(len)

for i,j in enumerate(x):
	x0 = x[i]
	a = y[i]
	ytmp = gaussian(x, x0, a=a, b=3.0)
	y2 = y2 + ytmp

plt.plot(x,y2)
plt.show()

