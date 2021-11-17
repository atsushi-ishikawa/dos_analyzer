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
#  python analyze.dos 1 Au_111 s p d
#
do_cohp = False
do_hilbert = False
check   = False
geometry_info = False

normalize_height  = True
relative_to_fermi = False

json = "surf_data.json"
db = connect(json)

argvs    = sys.argv
numpeaks = int(argvs[1])
system   = argvs[2]
doscar   = "DOSCAR_" + system
cohpcar  = "COHPCAR_" + system

print(" ----- %s ----- " % system)

sigma = 50

id  = db.get(system=system).id
obj = db[id]
atoms = db.get_atoms(id=id)  # surface + adsorbate

coord_num = get_adsorption_site(atoms=atoms, adsorbing_element="C")

orbitals = []
norbs = len(argvs) - 3  # number of orbitals
efermi = get_efermi_from_doscar(doscar)

# get orbital list
for i in range(0, norbs):
	orbitals.append(str(argvs[i+3]))

# finding natom
f = open(doscar, "r")
line1 = f.readline()
natom = int(line1.split()[0])
f.close()

vaspdos = VaspDos(doscar=doscar)
ene = vaspdos.energy
dos = vaspdos.dos

system  = obj.system
lattice = obj.lattice
data    = obj.data

# limit analysis on occupied part only
margin = 0.0
ene = list(filter(lambda x: x <= efermi+margin, ene))

for orbital in orbitals:
	# get pdos
	pdos = get_projected_dos(ene, natom, vaspdos, orbital)

	# smear
	pdos = smear_dos(ene, pdos, sigma=sigma)

	# get descriptor for pdos
	center, second = get_descriptors(range=range(0, natom))

	# get descriptor for surface
	center_surf, second_surf = get_descriptors(range=range(48, 64))

	# get descriptor for adsorption site
	center_site, second_site = get_descriptors(range=coord_ind)

	# analyze dos by fitting Gaussian
	peaks = findpeak(ene, pdos)

	width = 1.0*(1/sigma**0.5)  # guess for the Gaussian fit
	params = []
	for idx in peaks:
		params.append(ene[idx])
		params.append(pdos[idx])
		params.append(width)

	# Try gaussian fit. If fails, return blanked list
	try:
		params, rss, r2 = gaussian_fit(ene, pdos, params)
		peaks = sort_peaks(params, key="height")
		print("found %d peaks -- %s compoent R^2 = %5.3f" % (len(peaks), orbital, r2))
	except:
		params = []
		r2 = 0.0
		peaks = [(0, 0, 0) for i in range(numpeaks)]

	# discard if R^2 is too low
	if r2 < 0.90:  # 0.98:
		print("fitting failed: r2 (%5.3f) is too low ... quit" % r2)
		peaks = [(0, 0, 0) for i in range(numpeaks)]

	# sort peaks and limit to numpeaks
	peaks = peaks[0:numpeaks]
	tmp = []
	for peak in peaks:
		tmp.append(peak[0]); tmp.append(peak[1]); tmp.append(peak[2])

	peaks = sort_peaks(tmp, key="position")

	if relative_to_fermi:
		peaks = list(map(lambda x: (x[0]-efermi, x[1], x[2]), peaks))

	occ_peaks = list(filter(lambda x: x[0] < efermi+margin, peaks))
	vir_peaks = list(filter(lambda x: x[0] >= efermi+margin, peaks))

	# zero padding upto numpeaks
	occ_peaks = occ_peaks + [(0, 0, 0)]*(numpeaks - len(occ_peaks))
	vir_peaks = vir_peaks + [(0, 0, 0)]*(numpeaks - len(vir_peaks))

	# take upper edge by inverse Hilbert transform
	if do_hilbert:
		upper_edge, lower_edge = calc_edge(numpeaks, ene, pdos)
		upper_edge_surf, lower_edge_surf = calc_edge(numpeaks, ene, pdos_surf)
		upper_edge_site, lower_edge_site = calc_edge(numpeaks, ene, pdos_site)

	if do_cohp:
		cohp_pos_peak, cohp_pos_center, cohp_neg_peak, cohp_neg_center = cohp_analysis(cohpcar)

	# occupied and virtual
	position_occ = []; height_occ = []; width_occ = []; area_occ = []
	position_vir = []; height_vir = []; width_vir = []; area_vir = []

	for peak in occ_peaks:
		position_occ.append(peak[0])
		height_occ.append(peak[1])
		width_occ.append(peak[2])
		area_occ.append(peak[1]*peak[2]*np.sqrt(np.pi))
		# integral of Gauss function = height * sqrt(pi/a) where a = 1/width^2

	for peak in vir_peaks:
		position_vir.append(peak[0])
		height_vir.append(peak[1])
		width_vir.append(peak[2])
		area_vir.append(peak[1]*peak[2]*np.sqrt(np.pi))

	# if you want to check by eye
	if check:
		plot_dos(do_cohp, do_hilbert)

	# write to database
	data.update({orbital + "-" + "position_occ": position_occ})
	data.update({orbital + "-" + "height_occ": height_occ})
	data.update({orbital + "-" + "width_occ": width_occ})

	#upper_edge -= efermi
	#lower_edge -= efermi

if geometry_info:
	data = include_geometry_information(data)

if do_cohp:
	data.update({"cohp_pos_center": cohp_pos_center})
	data.update({"cohp_neg_center": cohp_neg_center})
	data.update({"cohp_pos_peak": cohp_pos_peak})
	data.update({"cohp_neg_peak": cohp_neg_peak})

db2 = connect("tmpout.json")
if relative_to_fermi:
	db2.write(atoms, system=system, lattice=lattice, data=data)
else:
	db2.write(atoms, system=system, lattice=lattice, data=data, efermi=efermi)

print("DONE for system =", system)

def include_geometry_information(data: dict):
	a, b, c, alpha, beta, gamma = atoms.get_cell_lengths_and_angles()
	surf_area = a * b * math.sin(math.radians(gamma))
	volume = atoms.get_volume()
	atom_num = atoms.numbers.sum()

	data.update({"surf_area": surf_area})
	data.update({"volume": volume})
	data.update({"atomic_numbers": atom_num})
	data.update({"atomic_number_site": coord_num})

	return data

def cohp_analysis(cohpcar):
	ene2, cohp = read_cohp(cohpcar)
	ene2 = np.array(ene2) + efermi

	cohp_pos = np.array(list(map(lambda x: x if x > 0.0 else 0.0, cohp)))
	cohp_neg = np.array(list(map(lambda x: x if x < 0.0 else 0.0, cohp)))

	cohp_pos = smear_dos(ene2, cohp_pos, sigma=sigma)
	cohp_neg = smear_dos(ene2, cohp_neg, sigma=sigma)

	# find peak for pos
	peaks  = findpeak(ene2, cohp_pos)
	maxind = peaks[np.argmax(cohp_pos[peaks])]  # highest peak
	cohp_pos_peak = ene2[maxind]

	# find peak for neg
	peaks = findpeak(ene2, cohp_neg)
	minind = peaks[np.argmin(cohp_neg[peaks])]  # highest peak
	cohp_neg_peak = ene2[minind]

	cohp_pos_center = get_distribution_moment(ene2, cohp_pos, order=1)
	cohp_neg_center = get_distribution_moment(ene2, cohp_neg, order=1)

	return cohp_pos_peak, cohp_pos_center, cohp_neg_peak, cohp_neg_center

def plot_dos(do_cohp=False, do_hilbert=False):
	"""

	Args:
		do_cohp:
		do_hilbert:

	Returns:

	"""
	import matplotlib.pylab as plt
	import seaborn as sb
	from scipy import fftpack

	fit = fit_func(ene, *params)

	if do_hilbert:
		ih = fftpack.ihilbert(pdos)

	sb.set(context='notebook', style='darkgrid', palette='deep',
		font='sans-serif', font_scale=1, color_codes=False, rc=None)
	plt.plot(ene, fit, label="fitted")
	plt.plot(ene, pdos, label="original")
	plt.vlines(position_occ, 0, np.max(pdos), linestyle="dashed", linewidth=0.5)
	# plt.plot(ene, ih, label="inverse Hilbert")
	plt.xlim([min(ene), max(ene) + margin])

	plt.legend()
	plt.show()

	if do_cohp:
		cohp_smeared = smear_dos(ene2, cohp, sigma=sigma)
		plt.plot(ene2, cohp_smeared, label="COHP")
		plt.legend()
		plt.show()

	return None

def get_descriptors(dos, range=None, orbital=None, normalize_height=True):
	"""
	get descriptors.

	Args:
		dos:
		range:
		orbital:
		normalize_height:
	Returns:

	"""
	tmp = np.zeros(len(ene))
	for i in range(0, natom):
		tmp += dos.site_dos(i, orbital)[:len(ene)]

	tmp = tmp[:len(ene)]
	tmp = tmp/np.max(tmpdos) if normalize_height else tmpdos

	center = get_distribution_moment(ene, tmp, order=1)
	second = get_distribution_moment(ene, tmp, order=2)
	third  = get_distribution_moment(ene, tmp, order=3)

	return center, second, third

def get_adsorption_site(atoms=None, adsorbing_element="C"):
	"""
	get coordinating atom.

	Args:
		atoms:
		adsorbing_element:

	Returns:

	"""
	coord_ind = []
	adsorbate_ind = atoms.get_chemical_symbols().index(adsorbing_element)
	li = [atoms.get_distance(adsorbate_ind, i) for i in range(len(atoms))]
	li = li[:adsorbate_ind]  # limit on surface
	li = list(map(lambda x: 1.0e3 if x < 1.0e-3 else x, li))

	coord_ind.append(np.argmin(li))
	coord_num = atoms[coord_ind].get_atomic_numbers()[0]  # atomic number of coordinating atom

	return coord_num

def get_projected_dos(ene, natom, vaspdos, orbital):
	"""
	Get projected DOS.

	Args:
		ene:
		natom:
		vaspdos: VaspDos object
		orbital:

	Returns:

	"""
	pdos = np.zeros(len(ene))
	for i in range(0, natom):
		pdos += vaspdos.site_dos(i, orbital)[:len(ene)]

	pdos = pdos[:len(ene)]
	return pdos

def normalize_height(dos):
	"""

	Args:
		dos:

	Returns:

	"""
	dos = dos / np.max(dos)
	return dos
