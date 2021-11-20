from ase import Atoms, Atom
from ase.calculators.vasp import VaspDos
import sys
import numpy as np
import math
from vasptools import *
from ase.db import connect
from vasptools import get_efermi_from_doscar
from ase.dft import get_distribution_moment
import matplotlib.pylab as plt

class VaspDosPlus:
	def __init__(self, doscar=None, numpeaks=1, system=None, orbitals=None):
		self._doscar = doscar
		self._numpeaks = numpeaks
		self._system = system
		self._jsonfile = None

		self.margin = 0.0

		self.vaspdos = VaspDos(doscar=doscar)
		self.efermi  = get_efermi_from_doscar(doscar)

		# limit analysis on occupied part only
		energy  = self.vaspdos.energy
		self.energy = list(filter(lambda x: x <= self.efermi + self.margin, energy))

		if orbitals is None:
			self.orbitals = {"s": 0, "p": 1, "d": 2}

		# finding natom
		with open(doscar, "r") as f:
			line1 = f.readline()
			self.natom = int(line1.split()[0])

		self._normalize = False
		self._do_hilbert = False
		self._do_cohp = False
		self._geometry_information = False

	def load_surface_data(self, json=None):
		self._jsonfile = json

	def get_descriptors(self):
		"""
		Get descriptors.

		Returns:
			dict

		Examples:
		python analyze.dos 1 Au_111 s p d
		"""
		check = False
		relative_to_fermi = False

		json = self._jsonfile
		db = connect(json)

		numpeaks = self._numpeaks
		system   = self._system
		doscar   = self._doscar
		if self._do_cohp:
			cohpcar  = self._system + "COHPCAR"

		print(" ----- %s ----- " % system)

		sigma = 50

		id  = db.get(system=system).id
		obj = db[id]
		atoms = db.get_atoms(id=id)  # surface + adsorbate

		# coord_num = get_adsorption_site(atoms=atoms, adsorbing_element="C")

		system  = obj.system
		lattice = obj.lattice

		descriptors = {}
		for orbital in self.orbitals.values():
			# get pdos
			pdos = self.get_projected_dos(self.vaspdos, orbital)

			# smear
			pdos = self.smear_dos(pdos, sigma=sigma)

			# get descriptor for pdos
			center, second = self.get_descriptor(atom_range=range(0, self.natom), orbital=orbital)

			# get descriptor for surface
			#center_surf, second_surf = get_descriptor(range=range(48, 64))

			# get descriptor for adsorption site
			#center_site, second_site = get_descriptor(range=coord_ind)

			# analyze dos by fitting Gaussian
			peaks = self.findpeak(pdos)

			width = 1.0*(1/sigma**0.5)  # guess for the Gaussian fit
			params = []
			for idx in peaks:
				params.append(self.energy[idx])
				params.append(pdos[idx])
				params.append(width)

			# Try gaussian fit. If fails, return blanked list
			try:
				params, rss, r2 = gaussian_fit(np.array(self.energy), pdos, params)
				peaks = sort_peaks(params, key="height")
				print("found %d peaks -- %s compoent R^2 = %5.3f" % (len(peaks), orbital, r2))
			except:
				r2 = 0.0
				peaks = [(0, 0, 0) for _ in range(numpeaks)]

			# discard if R^2 is too low
			if r2 < 0.90:  # 0.98:
				print("fitting failed: r2 (%5.3f) is too low ... quit" % r2)
				peaks = [(0, 0, 0) for _ in range(numpeaks)]

			# sort peaks and limit to numpeaks
			peaks = peaks[0:numpeaks]
			tmp = []
			for peak in peaks:
				tmp.append(peak[0]); tmp.append(peak[1]); tmp.append(peak[2])

			peaks = sort_peaks(tmp, key="position")

			if relative_to_fermi:
				peaks = list(map(lambda x: (x[0] - self.efermi, x[1], x[2]), peaks))

			occ_peaks = list(filter(lambda x: x[0] < self.efermi + self.margin, peaks))
			vir_peaks = list(filter(lambda x: x[0] >= self.efermi + self.margin, peaks))

			# zero padding upto numpeaks
			occ_peaks = occ_peaks + [(0, 0, 0)]*(numpeaks - len(occ_peaks))
			vir_peaks = vir_peaks + [(0, 0, 0)]*(numpeaks - len(vir_peaks))

			# take upper edge by inverse Hilbert transform
			if self._do_hilbert:
				upper_edge, lower_edge = calc_edge(numpeaks, pdos)

			if self._do_cohp:
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
				plot_dos()

			# write to database
			orbitalname = get_key_from_value(self.orbitals, orbital)
			descriptors.update({orbitalname + "-" + "position_occ": position_occ})
			descriptors.update({orbitalname + "-" + "height_occ": height_occ})
			descriptors.update({orbitalname + "-" + "width_occ": width_occ})
			descriptors.update({orbitalname + "-" + "center": center})
			descriptors.update({orbitalname + "-" + "second": second})

		#upper_edge -= efermi
			#lower_edge -= efermi

		if self._geometry_information:
			descriptors = self.include_geometry_information(atoms, descriptors)

		if self._do_cohp:
			descriptors.update({"cohp_pos_center": cohp_pos_center})
			descriptors.update({"cohp_neg_center": cohp_neg_center})
			descriptors.update({"cohp_pos_peak": cohp_pos_peak})
			descriptors.update({"cohp_neg_peak": cohp_neg_peak})

		if self._do_hilbert:
			descriptors.update({"upper_edge": upper_edge})
			descriptors.update({"lower_edge": lower_edge})

		print("DONE for system =", system)
		return descriptors

	def include_geometry_information(self, atoms, descriptors: dict):
		a, b, c, alpha, beta, gamma = atoms.get_cell_lengths_and_angles()
		surf_area = a * b * math.sin(math.radians(gamma))
		volume = atoms.get_volume()
		atom_num = atoms.numbers.sum()

		descriptors.update({"surf_area": surf_area})
		descriptors.update({"volume": volume})
		descriptors.update({"atomic_numbers": atom_num})
		descriptors.update({"atomic_number_site": coord_num})

		return data

	def cohp_analysis(self, cohpcar):
		ene2, cohp = read_cohp(cohpcar)
		ene2 = np.array(ene2) + self.efermi

		cohp_pos = np.array(list(map(lambda x: x if x > 0.0 else 0.0, cohp)))
		cohp_neg = np.array(list(map(lambda x: x if x < 0.0 else 0.0, cohp)))

		cohp_pos = smear_dos(cohp_pos, sigma=sigma)
		cohp_neg = smear_dos(cohp_neg, sigma=sigma)

		# find peak for pos
		peaks  = findpeak(cohp_pos)
		maxind = peaks[np.argmax(cohp_pos[peaks])]  # highest peak
		cohp_pos_peak = ene2[maxind]

		# find peak for neg
		peaks = findpeak(cohp_neg)
		minind = peaks[np.argmin(cohp_neg[peaks])]  # highest peak
		cohp_neg_peak = ene2[minind]

		cohp_pos_center = get_distribution_moment(cohp_pos, order=1)
		cohp_neg_center = get_distribution_moment(cohp_neg, order=1)

		return cohp_pos_peak, cohp_pos_center, cohp_neg_peak, cohp_neg_center

	def plot_dos(self):
		"""
		Plot dos.
		"""
		import seaborn as sb
		from scipy import fftpack

		fit = fit_func(self.energy, *params)

		if self._do_hilbert:
			ih = fftpack.ihilbert(pdos)

		sb.set(context='notebook', style='darkgrid', palette='deep',
			font='sans-serif', font_scale=1, color_codes=False, rc=None)
		plt.plot(self.energy, fit, label="fitted")
		plt.plot(self.energy, pdos, label="original")
		plt.vlines(position_occ, 0, np.max(pdos), linestyle="dashed", linewidth=0.5)
		if self._do_hilbert:
			plt.plot(self.energy, ih, label="inverse Hilbert")
		plt.xlim([min(self.energy), max(self.energy) + self.margin])

		plt.legend()
		plt.show()

		if self._do_cohp:
			cohp_smeared = smear_dos(ene2, cohp, sigma=sigma)
			plt.plot(ene2, cohp_smeared, label="COHP")
			plt.legend()
			plt.show()

		return None

	def get_descriptor(self, atom_range=None, orbital=None):
		"""
		get descriptors.

		Args:
			atom_range:
			orbital:
		Returns:

		"""
		tmp = np.zeros(len(self.energy))
		for i in atom_range:
			tmp += self.vaspdos.site_dos(i, orbital)[:len(self.energy)]

		tmp = tmp[:len(self.energy)]
		tmp = self.normalize_height(tmp) if self._normalize else tmp

		center = get_distribution_moment(self.energy, tmp, order=1)
		second = get_distribution_moment(self.energy, tmp, order=2)

		return center, second

	def get_adsorption_site(self, atoms=None, adsorbing_element="C"):
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

	def get_projected_dos(self, vaspdos, orbital):
		"""
		Get projected DOS.

		Args:
			vaspdos: VaspDos object
			orbital:

		Returns:

		"""
		pdos = np.zeros(len(self.energy))
		for i in range(0, self.natom):
			pdos += vaspdos.site_dos(i, orbital)[:len(self.energy)]

		pdos = pdos[:len(self.energy)]
		return pdos

	def normalize_height(self, dos):
		"""

		Args:
			dos:

		Returns:

		"""
		dos = dos / np.max(dos)
		return dos

	def gaussian(self, x, x0, a, b):
		"""
		y = a*exp(-b*(x-x0)**2)
		"""
		import numpy as np
		x = np.array(x)
		y = np.exp(-b*(x-x0)**2)
		y = a*y
		return y

	def smear_dos(self, dos, sigma=5.0):
		"""
		get smeared dos.

		Args:
			dos:
			sigma:

		Returns:

		"""
		x = self.energy
		y = dos
		smeared = np.zeros(len(x))

		for i, j in enumerate(x):
			x0 = x[i]
			a = y[i]
			smeared += self.gaussian(x, x0, a=a, b=sigma)

		return smeared

	def findpeak(self, y):
		import peakutils

		indexes = peakutils.indexes(y, thres=0.1, min_dist=1)
		return indexes

def gaussian_fit(x, y, guess):
	from scipy.optimize import curve_fit

	def fit_func(x, *params):
		y = np.zeros_like(x)
		for i in range(0, len(params), 3):
			ctr = params[i]  # position
			amp = params[i + 1]  # height
			wid = params[i + 2]  # width (smaller is sharper)
			y = y + amp * np.exp(-((x - ctr) / wid) ** 2)
		return y

	#popt, pcov = curve_fit(fit_func, x, y, p0=guess, method="trf", ftol=1.0e-5, xtol=1.0e-5)
	#popt, pcov = curve_fit(fit_func, x, y, p0=guess, method="trf", ftol=1.0e-6, xtol=1.0e-6)
	popt, pcov = curve_fit(fit_func, x, y, p0=guess, method="lm", ftol=1.0e-5, xtol=1.0e-5)  # good
	#popt, pcov = curve_fit(fit_func, x, y, p0=guess, method="lm", ftol=1.0e-6, xtol=1.0e-6)

	fit = fit_func(x, *popt)
	residual = y - fit
	rss = np.sum(residual**2)  # residual sum of squares
	tss = np.sum((y-np.mean(y))**2)  # total sum of squares
	r2 = 1 - (rss / tss)

	return popt, rss, r2

def sort_peaks(peaks, key="height"):
	"""
	assuming peaks are stored in [position, height, width,  position, height, width,...]
	"""
	dtype = [("position", float), ("height", float), ("width", float)]

	newpeaks = np.array([], dtype=dtype)

	for i in range(0, len(peaks), 3):
		peak = np.array((peaks[i], peaks[i+1], peaks[i+2]), dtype=dtype)
		newpeaks = np.append(newpeaks, peak)

	newpeaks = np.sort(newpeaks, order=key)
	newpeaks = newpeaks[::-1]  # sort in descending order
	return newpeaks

def get_key_from_value(d, val):
	keys = [k for k, v in d.items() if v == val]
	if keys:
		return keys[0]
