from ase import Atoms, Atom
from ase.calculators.vasp import VaspDos
import numpy as np
from vasptools import *
from ase.db import connect

json  = "tmp2.json"
db = connect(json)
system = "Pd111"
doscar = "DOSCAR_" + system
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

#id = db.get(system=system).id
#print id
#db.update(id, peak0=peaks[0])

id  = db.get(system=system).id
obj = db[id]

system   = obj.system # already know!
lattice  = obj.lattice
olddata = obj.data

pos = [peaks[0][0],peaks[1][0]]
tmpdata  = {"pos": pos}
newdata = olddata + tmpdata

atoms = db.get_atoms(id=id)

db2 = connect("tmpout.json")
db2.write(atoms,system=system,lattice=lattice,data=newdata)

