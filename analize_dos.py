from ase import Atoms, Atom
from ase.calculators.vasp import VaspDos
import numpy as np
from vasptools import *
from ase.db import connect

json = "surf.json"
db = connect(json)
doscar = "DOSCAR_Pd111"
sigma = 100.0

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

width = 0.05
params = []
for idx in peaks:
 	params.append(ene[idx])
 	params.append(ddos[idx])
	params.append(width)
	
params = gaussian_fit(ene, ddos, params)
peaks = sort_peaks_by_height(params)
print peaks
#
# checking by eye
#
#import matplotlib.pylab as plt
#from vasptools import fit_func
#fit = fit_func(ene,*params)
#plt.plot(ene,fit)
#plt.plot(ene,ddos)
#plt.show()

id = db.get(element=element,face=face_str).id
db.update(id,dos=...)

