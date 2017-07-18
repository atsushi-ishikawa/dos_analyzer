def gaussian(x,x0,a,b):
	"""
	  y = a*exp(-b*(x-x0)**2)
	"""
	import numpy as np
	x = np.array(x)
	y = np.exp( -b*(x-x0)**2 )
	y = a*y
	return y

def smear_dos(dos, sigma=5.0):
	"""
	  get smeared dos
	"""
	import numpy as np

	x = dos.energy
	y = dos.dos

	len = len(x)
	y2  = np.zeros(len)

	for i,j in enumerate(x):
		x0 = x[i]
		a = y[i]
		ytmp = gaussian(x, x0, a=a, b=5.0)
		y2 = y2 + ytmp

	return y2

