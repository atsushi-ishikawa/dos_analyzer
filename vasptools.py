def lattice_info_guess(bulk):
	from ase import Atoms
	#
	# dict for already-known cases
	#
	dict = {
		"Rh" : {"lattice" : "fcc" , "a0" : 3.803},
		"Pd" : {"lattice" : "fcc" , "a0" : 3.891},
		"Ir" : {"lattice" : "fcc" , "a0" : 3.839},
		"Pt" : {"lattice" : "fcc" , "a0" : 3.924},
		"Cu" : {"lattice" : "fcc" , "a0" : 3.615},
		"Ag" : {"lattice" : "fcc" , "a0" : 4.085},
		"Au" : {"lattice" : "fcc" , "a0" : 4.078}
	}
	element = bulk.get_chemical_symbols()[0]
	if element in dict:
		print("I know this.")
		lattice = dict[element]["lattice"]
		a0 = dict[element]["a0"]
	else:
		print("I don't know this.")
		lattice = "fcc"
		a0 = 4.0

	return lattice, a0

def make_bulk(element1, element2=None, comp1=100, lattice="fcc", a0=4.0, repeat=1):
	from ase import Atoms
	from ase.build import bulk
	import numpy as np
	from ase.visualize import view
	#
	# to fractional
	#
	if element2 is not None:
		comp2 = 100 - comp1
		comp1 = comp1/100.0
		comp2 = comp2/100.0
	#
	# form bulk first
	#
	bulk = bulk(element1, lattice, a=a0, cubic=True)
	bulk = bulk.repeat( repeat ) # for consistency with alloys

	if element2 is not None:
		#
		# make alloy if needed
		#
		natom_tot = len(bulk.get_atomic_numbers())
		natom2    = int(comp2 * natom_tot)
		natom1    = natom_tot - natom2
		#
		# replace atoms in bulk randomly
		#
		list = bulk.get_chemical_symbols()
		#
		#np.random.seed(1) # set random seed for reproducability --> gives same structures for all alloys
		np.random.seed() 

		replace_list = np.random.choice(len(list), size=natom2, replace=False)
		for i in replace_list:
			list[i] = element2

		bulk.set_chemical_symbols(list)

	return bulk

def optimize_lattice_constant(bulk, lattice="fcc", a0=4.0, xc="PBEsol", encut=400, ediff=1.0e-5, ediffg=1.0e-6, npar=1, nsim=1):
	""" 
	function to do bulk optimization
	"""
	from ase import Atoms
	from ase.calculators.vasp import Vasp
	import os, shutil
	#
	# directry things
	#
	cudir   = os.getcwd()
	workdir = os.path.join(cudir, "lattice_opt")
	os.makedirs(workdir)
	os.chdir(workdir)
	os.system("mkdir work_bulk")
	#
	# compuational condition for Vasp
	#
	prec   = "normal"
	potim  = 0.1
	nelmin = 5
	kpts   = [3,3,3]
	gamma  = True
	nsw    = 200
	isif   = 6 # or 6---freezing ions
	isym   = -1

	xc = xc.lower()
	if xc == "pbe" or xc == "pbesol" or xc == "rpbe":
		pp = "pbe"
	elif xc == "pw91":
		pp = "pw91"
	elif xc == "lda":
		pp = "lda"
	else:
		print("xc error")

	calc = Vasp(prec=prec, xc=xc, pp=pp, ispin=1,
				ismear=1, sigma=0.2, isif=isif, nelmin=nelmin, encut=encut,
				ibrion=2, nsw=nsw, potim=potim, ediff=ediff, ediffg=ediffg,
				kpts=kpts, gamma=gamma, isym=isym, npar=npar, nsim=nsim, lreal=False )

	bulk.set_calculator(calc)
	bulk.get_potential_energy()

	cell = bulk.cell # optimized cell size

	os.chdir(cudir)
	shutil.rmtree(workdir)

def gaussian(x,x0,a,b):
	"""
	  y = a*exp(-b*(x-x0)**2)
	"""
	import numpy as np
	x = np.array(x)
	y = np.exp( -b*(x-x0)**2 )
	y = a*y
	return y

def smear_dos(energy, dos, sigma=5.0):
	"""
	  get smeared dos
	"""
	import numpy as np

	x = energy
	y = dos

	length = len(x)
	y2  = np.zeros(length)

	for i,j in enumerate(x):
		x0 = x[i]
		a = y[i]
		ytmp = gaussian(x, x0, a=a, b=sigma)
		y2 = y2 + ytmp

	return y2

def sort_atoms_by_z(atoms):
	from ase import Atoms, Atom
	import numpy as np
	#
	# keep information for original Atoms
	#
	tags = atoms.get_tags()
	pbc  = atoms.get_pbc()
	cell = atoms.get_cell()

	dtype = [("idx",int), ("z",float)]
	zlist = np.array([], dtype=dtype)

	for idx, atom in enumerate(atoms):
		tmp = np.array([(idx,atom.z)],dtype=dtype)
		zlist = np.append(zlist, tmp)

	zlist = np.sort(zlist, order="z")

	newatoms = Atoms()

	for i in zlist:
		idx = i[0]
		newatoms.append(atoms[idx])
	#
	# restore
	#
	newatoms.set_tags(tags)
	newatoms.set_pbc(pbc)
	newatoms.set_cell(cell)

	return newatoms

def findpeak(x, y):
	import numpy as np
	import peakutils 

	indexes = peakutils.indexes(y, thres=0.1, min_dist=1)
	return indexes

def fit_func(x, *params):
	import numpy as np
	y = np.zeros_like(x)
	for i in range(0, len(params), 3):
		ctr = params[i]     # position
		amp = params[i+1]   # height
		wid = params[i+2]   # width (smaller is sharper)
		y   = y + amp * np.exp( -((x-ctr)/wid)**2 )
	return y

def gaussian_fit(x, y, guess):
	from scipy.optimize import curve_fit
	import numpy as np

	popt, pcov = curve_fit(fit_func, x, y, p0=guess)

	fit = fit_func(x,*popt)
	residual = y - fit
	rss = np.sum(residual**2) # residual sum of squares
	tss = np.sum( (y-np.mean(y))**2 ) # total sum of squares
	r2 = 1 - (rss / tss)

	return popt,rss,r2

def sort_peaks(peaks, key="height"):
	import numpy as np
	"""
	  assuming peaks are stored in [position, height, width,  position, height, width,...]
	"""
	dtype = [ ("position",float), ("height",float), ("width",float) ]

	newpeaks = np.array([], dtype=dtype)

	for i in range(0, len(peaks), 3):
		peak = np.array(( peaks[i],peaks[i+1],peaks[i+2] ), dtype=dtype)
		newpeaks = np.append(newpeaks, peak)

	newpeaks = np.sort(newpeaks, order=key)
	newpeaks = newpeaks[::-1] # sort in descending order
	return newpeaks

def get_efermi_from_doscar(DOSCAR):
	import linecache
	line   = linecache.getline(DOSCAR,6)
	line   = line.strip()
	line   = line.split()
	efermi = float(line[3])
	return efermi

def calc_edge(numpeaks, ene, pdos):
	import numpy as np
	from scipy import fftpack
	#
	# take upper edge by inverse Hilbert transform
	#
	ih = fftpack.ihilbert(pdos)
	ihpeaks = findpeak(ene, abs(ih))

	upper_edge = np.zeros(numpeaks)
	lower_edge = np.zeros(numpeaks)

	for peak in ihpeaks:
		if ih[peak] > 0.0 and ih[peak] > 0.8*max(ih): # just choose large peak, positive part
			upper_edge = np.insert(upper_edge, 0, ene[peak])
		elif ih[peak] <= 0.0 and ih[peak] < 0.8*min(ih):
			lower_edge = np.insert(lower_edge, 0, ene[peak])

	upper_edge = upper_edge[0:numpeaks]
	lower_edge = lower_edge[0:numpeaks]
	#upper_edge = upper_edge[::-1]
	#lower_edge = lower_edge[::-1]

	return upper_edge, lower_edge

def make_lobsterin():
	"""
	make input file for lobster (lobsterin) for alloy systems
	"""
	f = open("lobsterin", "w")
	str = ["COHPstartEnergy -20", "COHPendEnergy 20", "basisSet pbeVaspFit2015", "includeOrbitals spd", "cohpbetween atom 1 and atom 2"]
	str = " \n".join(str)
	f.write(str)
	f.close()

def read_cohp(cohpfile="COHPCAR.lobster"):
	f = open(cohpfile, "r")
	f.readline() # skip first
	ndos = int(f.readline().split()[2])

	[f.readline() for i in range(2)]

	ene  = []
	cohp = []

	for i in range(ndos):
		tmp = f.readline().split()
		ene.append( float(tmp[0]))
		cohp.append(float(tmp[1]))
	
	f.close()
	return ene,cohp

