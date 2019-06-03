from ase import Atoms, Atom
from ase.calculators.vasp import VaspDos
import sys
import numpy as np
from vasptools import *
from ase.db import connect
from vasptools import fit_func
from vasptools import get_efermi_from_doscar
#
# Usage:
#  python analyze.dos Cu_111 s p d
#
json  = "surf_data_alloy.json"
db = connect(json)

argvs  = sys.argv
doscar = argvs[1]
system = doscar.split("_")[1] + "_" + doscar.split("_")[2]

sigma  = 2.0 # smaller the broader. 0.5 gives broaden peak -- for singple peak use
check  = False
numpeaks = 1

orbitals = []
norbs = len(argvs) - 2 # number of orbitals

efermi = get_efermi_from_doscar(doscar)
#
# get orbital list
#
for i in range(0, norbs):
	orbitals.append(str(argvs[i+2]))
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

margin = 1.0
ene2 = list(filter(lambda x : x <= efermi+margin, ene))
old_last = len(ene)
new_last = len(ene2)
ene = ene2

id  = db.get(system=system).id
obj = db[id]

system  = obj.system # already know!
lattice = obj.lattice
data    = obj.data

for orbital in orbitals:
	tmppdos = np.zeros(old_last)
	pdos = np.zeros(new_last)
	for i in range(0, natom):
		tmppdos = tmppdos + dos.site_dos(i, orbital)
	pdos = tmppdos[:new_last]

	pdos = smear_dos(ene, pdos, sigma=sigma)
	#pdos = pdos/max(pdos) # relative
	peaks = findpeak(ene, pdos)
	if check:
		for i in peaks:
			print(ene[i])

	width = 0.10 # guess for the Gaussian fit
	params = []
	for idx in peaks:
		params.append(ene[idx])
		params.append(pdos[idx])
		params.append(width)
	#
	# Try gaussian fit. If fails, return blanked list
	#
	try:
		params = gaussian_fit(ene, pdos, params)
		peaks  = sort_peaks(params,key="position") # sort by pos. would be better than height
	except:
		params = []
		peaks  = [(0,0,0) for i in range(numpeaks)]
	#
	# if you want to check by eye
	#
	if check:
		import matplotlib.pylab as plt
		import seaborn as sb

		print("height  : ", end="")
		for i in range(numpeaks):
			print("%12.6f" % peaks[i][1], end=" ")
		print("")
		print("position: ", end="")
		for i in range(numpeaks):
			print("%12.6f" % peaks[i][0], end=" ");
		print("")

		fit = fit_func(ene,*params)
		sb.set(context='notebook', style='darkgrid', palette='deep', font='sans-serif', font_scale=1, color_codes=False, rc=None)
		plt.plot(ene,fit,  label="fitted")
		plt.plot(ene,pdos, label="original")
		plt.legend()
		plt.show()
	#
	# adding to database
	#
	position = []
	height   = []
	width    = []

	for i in range(numpeaks):
		position.append(peaks[i][0])
		height.append(peaks[i][1])
		width.append(peaks[i][2])
	#
	# sort by position : only useful in numpeaks = 2
	#
	if numpeaks==2 and position[0] > position[1]:
#	if numpeaks==2 and height[1] > height[0]:
	 	tmp1 = position[0];	position[0] = position[1];  position[1] = tmp1
 		tmp2 = height[0];	height[0]   = height[1];	height[1]   = tmp2
 		tmp3 = width[0];	width[0]    = width[1];		width[1]    = tmp3

	data.update({ orbital + "-dos " + "position" : position})
	data.update({ orbital + "-dos " + "height"   : height})
	data.update({ orbital + "-dos " + "width"    : width})

	atoms = db.get_atoms(id=id)

db2 = connect("tmpout.json")
db2.write(atoms, system=system, lattice=lattice, data=data, efermi=efermi)

print("DONE for system =", system)

