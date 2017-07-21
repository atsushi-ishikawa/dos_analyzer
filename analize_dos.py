from ase import Atoms, Atom
from ase.calculators.vasp import VaspDos
import sys
import numpy as np
from vasptools import *
from ase.db import connect

json  = "surf_data.json"
db = connect(json)
#system = "Pd111"
argvs = sys.argv
system  = argvs[1]
orbital = argvs[2]
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

pdos = np.zeros(len(tdos))

for i in range(0,natom):
	pdos = pdos + dos.site_dos(i, orbital)

pdos = smear_dos(ene, pdos, sigma=sigma)

peaks = findpeak(ene, pdos)

width = 0.05
params = []
for idx in peaks:
 	params.append(ene[idx])
 	params.append(pdos[idx])
	params.append(width)
	
params = gaussian_fit(ene, pdos, params)
peaks = sort_peaks_by_height(params)
#
# checking by eye
#
#import matplotlib.pylab as plt
#from vasptools import fit_func
#fit = fit_func(ene,*params)
#plt.plot(ene,fit)
#plt.plot(ene,pdos)
#plt.show()

#
# adding to database
#
id  = db.get(system=system).id
obj = db[id]

system   = obj.system # already know!
lattice  = obj.lattice
data     = obj.data

numpeaks = 5
position = []
height   = []
width    = []
for i in range(numpeaks):
	position.append(peaks[i][0])
	height.append(peaks[i][1])
	width.append(peaks[i][2])

data.update({ orbital + "-dos " + "position" : position})
data.update({ orbital + "-dos " + "height"   : height})
data.update({ orbital + "-dos " + "width"    : width})

atoms = db.get_atoms(id=id)

db2 = connect("tmpout.json")
db2.write(atoms,system=system,lattice=lattice,data=data)

