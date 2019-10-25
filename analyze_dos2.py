from ase import Atoms, Atom
from ase.calculators.vasp import VaspDos
import sys
import numpy as np
import math
from vasptools import *
from ase.db import connect
from vasptools import fit_func
from vasptools import get_efermi_from_doscar
from ase.dft import get_distribution_moment

import matplotlib.pyplot as plt
#
# Usage:
#  python analyze.dos 1 DOSCAR_Au_111 s p d
#
json  = "surf_data.json"
db = connect(json)

argvs  = sys.argv
numpeaks = int(argvs[1])
doscar = argvs[2]
system = doscar.split("_")[1] + "_" + doscar.split("_")[2]

if numpeaks==1:
	sigma = 3.0 # 3.0 -- d-dominant but not good learning curve # 1.0,2.0 -- d not dominant
elif numpeaks==2:
	sigma = 5.0 # 2.0 # 1.0
else:
	sigma = 5.0 # 6.0 # 3.0 # 4.0

check = False

id  = db.get(system=system).id
obj = db[id]
atoms = db.get_atoms(id=id) # surface + adsorbate
#
# get coordinating atom
#
adsorbate = "C"
coord_ind = []
adsorbate_ind = atoms.get_chemical_symbols().index(adsorbate)
li = [atoms.get_distance(adsorbate_ind, i) for i in range(len(atoms)) ]
li = li[:adsorbate_ind] # limit on surface
li = list(map(lambda x : 1.0e3 if x<1.0e-3 else x, li))
coord_ind.append(np.argmin(li))

orbitals = []
norbs = len(argvs) - 3 # number of orbitals 
efermi = get_efermi_from_doscar(doscar)
#
# get orbital list
#
for i in range(0, norbs):
	orbitals.append(str(argvs[i+3]))
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

system  = obj.system # already know!
lattice = obj.lattice
data    = obj.data

for orbital in orbitals:
	#
	# get pdos for total system
	#
	tmppdos = np.zeros(len(ene))
	pdos    = np.zeros(len(ene))
	for i in range(0, natom):
		tmppdos = tmppdos + dos.site_dos(i,orbital)[:len(ene)]

	pdos   = tmppdos[:len(ene)]
	pdos   = smear_dos(ene, pdos, sigma=sigma)
	pdos   = pdos/np.max(pdos)
	peaks  = findpeak(ene, pdos)
	integ  =  get_distribution_moment(ene, tmppdos, order=0)
	center =  get_distribution_moment(ene, tmppdos, order=1)
	#
	# get pdos for surface layer
	#
	tmppdos   = np.zeros(len(ene))
	pdos_surf = np.zeros(len(ene))
	for i in range(48,64):
		tmppdos = tmppdos + dos.site_dos(i, orbital)[:len(ene)]
	pdos_surf   = tmppdos[:len(ene)]
	pdos_surf   = smear_dos(ene, pdos_surf, sigma=sigma)
	integ_surf  = get_distribution_moment(ene, tmppdos, order=0)
	center_surf = get_distribution_moment(ene, tmppdos, order=1)
	#
	# get pdos for adsorbating atom
	#
	tmppdos   = np.zeros(len(ene))
	pdos_site = np.zeros(len(ene))
	for i in coord_ind:
		tmppdos = tmppdos + dos.site_dos(i, orbital)[:len(ene)]
	pdos_site   = tmppdos[:len(ene)]
	pdos_site   = smear_dos(ene, pdos_site, sigma=sigma)
	integ_site  = get_distribution_moment(ene, tmppdos, order=0)
	center_site = get_distribution_moment(ene, tmppdos, order=1)

	width = 0.05 # guess for the Gaussian fit
	params = []
	for idx in peaks:
		params.append(ene[idx])
		params.append(pdos[idx])
		params.append(width)
	#
	# Try gaussian fit. If fails, return blanked list
	#
	try:
		params,rss,r2 = gaussian_fit(ene, pdos, params)
		peaks = sort_peaks(params, key="position")
		#peaks = list(map(lambda x : (x[0]-efermi,x[1],x[2]), peaks)) # relative from fermi
		print("R^2 = %5.3f" % r2)
	except:
		params = []
		r2 = 0.0
		peaks  = [(0,0,0) for i in range(numpeaks)]

	# quit if R^2 is too low
	if r2 < 0.0:
		print("fitting failed ... quit")
		peaks  = [(0,0,0) for i in range(numpeaks)]

	margin = 0.5
	occ_peaks = list(filter(lambda x : x[0] <  efermi+margin, peaks))
	vir_peaks = list(filter(lambda x : x[0] >= efermi+margin, peaks))

	# zero padding upto numpeaks
	occ_peaks = occ_peaks + [(0,0,0)]*(numpeaks - len(occ_peaks))
	vir_peaks = vir_peaks + [(0,0,0)]*(numpeaks - len(vir_peaks))
	#
	# take upper edge by inverse Hilbert transform
	#
	upper_edge, lower_edge = calc_edge(numpeaks, ene, pdos)
	upper_edge_surf, lower_edge_surf = calc_edge(numpeaks, ene, pdos_site)
	upper_edge_site, lower_edge_site = calc_edge(numpeaks, ene, pdos_site)
	#edge = np.hstack( (upper_edge, lower_edge) )
	#
	# if you want to check by eye
	#
	if check:
		import matplotlib.pylab as plt
		import seaborn as sb

		print("upper edge:{}".format(upper_edge))
		print("lower edge:{}".format(lower_edge))

		fit = fit_func(ene,*params)

		sb.set(context='notebook', style='darkgrid', palette='deep', font='sans-serif', font_scale=1, color_codes=False, rc=None)
		plt.plot(ene, fit,  label="fitted")
		plt.plot(ene, pdos, label="original")
		plt.plot(ene, ih, label="inverse Hilbert")
		plt.legend()
		plt.show()
	#
	# adding to database
	#
	position_occ = [] ; height_occ = [] ; width_occ = [] ; area_occ = []
	position_vir = [] ; height_vir = [] ; width_vir = [] ; area_vir = []

	for peak in occ_peaks:
		position_occ.append(peak[0])
		height_occ.append(peak[1])
		width_occ.append(peak[2])
		area_occ.append(peak[1]*peak[2]*np.sqrt(np.pi)) # integral of Gauss function = height * sqrt(pi/a) where a = 1/width^2
	for peak in vir_peaks:
		position_vir.append(peak[0])
		height_vir.append(peak[1])
		width_vir.append(peak[2])
		area_vir.append(peak[1]*peak[2]*np.sqrt(np.pi))

	data.update({ orbital + "-" + "position_occ" : position_occ})
	data.update({ orbital + "-" + "height_occ"   : height_occ})
	data.update({ orbital + "-" + "width_occ"    : width_occ})
	#data.update({ orbital + "-" + "area_occ"     : area_occ})

	data.update({ orbital + "-" + "position_vir" : position_vir})
	data.update({ orbital + "-" + "height_vir"   : height_vir})
	data.update({ orbital + "-" + "width_vir"    : width_vir})
	#data.update({ orbital + "-" + "area_vir"     : area_vir})

	upper_edge -= efermi
	lower_edge -= efermi

	#data.update({ orbital + "-" + "upper_edge" : upper_edge})
	#data.update({ orbital + "-" + "lower_edge" : lower_edge})
	#data.update({ orbital + "-" + "upper_edge_site" : upper_edge_site})
	#data.update({ orbital + "-" + "upper_edge_surf" : upper_edge_surf})

	#data.update({ orbital + "-" + "center"      : center})
	#data.update({ orbital + "-" + "center_surf" : center_surf})
	#data.update({ orbital + "-" + "center_site" : center_site})
	#data.update({ orbital + "-" + "integ"      : integ})
	#data.update({ orbital + "-" + "integ_site" : integ_site})
	#data.update({ orbital + "-" + "integ_surf" : integ_surf})

#a,b,c,alpha,beta,gamma = atoms.get_cell_lengths_and_angles()
#surf_area = a*b*math.sin(math.radians(gamma))
#volume = atoms.get_volume()
#data.update({"surf_area" : surf_area})
#data.update({"volume" : volume})
#atom_num = atoms.numbers.sum()
#data.update({"atomic_numbers" : atom_num})

db2 = connect("tmpout.json")
db2.write(atoms, system=system, lattice=lattice, data=data, efermi=efermi)
#db2.write(atoms, system=system, lattice=lattice, data=data)

print("DONE for system =", system)

