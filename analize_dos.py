from ase import Atoms, Atom
from ase.calculators.vasp import VaspDos
import numpy as np
import matplotlib.pylab as plt
from peakutils.plot import plot as pplot

from vasptools import smear_dos

from tmp import *

doscar = "DOSCAR_Pd111"
sigma = 10.0

#
# finding natom
#
f = open(doscar, "r")
line1 = f.readline()
natom = int( line1.split()[0] )
f.close()

dos = VaspDos(doscar=doscar)

ene  = dos.energy
tdos = dos.dos

sdos = np.zeros(len(tdos))
pdos = np.zeros(len(tdos))
ddos = np.zeros(len(tdos))

for i in range(0,natom):
	sdos = sdos + dos.site_dos(i, "s")
	pdos = pdos + dos.site_dos(i, "p")
	ddos = ddos + dos.site_dos(i, "d")

tdos = smear_dos(ene, tdos, sigma=sigma)
sdos = smear_dos(ene, sdos, sigma=sigma)
pdos = smear_dos(ene, pdos, sigma=sigma)
ddos = smear_dos(ene, ddos, sigma=sigma)

peaks = findpeak(ene, ddos)

width = 0.1
params = []
for idx in peaks:
 	params.append(ene[idx])
 	params.append(ddos[idx])
	params.append(width)
	
print params
params = gaussian_fit(ene, ddos, params)
print params
#plt.plot(ene,sdos)
#plt.plot(ene,pdos)
pplot(ene,ddos,peaks)
plt.show()

