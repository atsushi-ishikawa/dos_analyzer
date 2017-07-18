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
		print "I know this."
		lattice = dict[element]["lattice"]
		a0 = dict[element]["a0"]
	else:
		print "I don't know this."
		lattice = "fcc"
		a0 = 4.0

	return lattice,a0

def make_bulk(element1, element2=None, comp1=100,lattice="fcc",a0=4.0,repeat=2):
	from ase import Atoms
	from ase.build import bulk
	import numpy as np
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

	if element2 is not None:
		#
		# make alloy if needed
		#
		bulk = bulk.repeat( (repeat,repeat,repeat) )
		natom_tot = len(bulk.get_atomic_numbers())
		natom2    = int(comp2 * natom_tot)
		natom1    = natom_tot - natom2
		#
		# replace atoms in bulk randomly
		#
		list = bulk.get_chemical_symbols()
		#
		replace_list = np.random.choice(len(list), size=natom2)
		for i in replace_list:
			list[i] = element2

		bulk.set_chemical_symbols(list)

	return bulk

def get_optimized_lattice_constant(bulk, lattice="fcc",a0=4.0, xc="PBEsol"):
	""" function to return optimized bulk constant
	"""
	from ase import Atoms
#	from ase.calculators.vasp import Vasp
 	from ase.calculators.emt import EMT
	#
	# compuational condition for Vasp
	#
	prec   = "normal"
	potim  = 0.1
	ediff  = 1.0e-6
	ediffg = -0.001
	kpts = [10, 10, 10]

	xc = xc.lower()
	if xc == "pbe" or xc == "pbesol" or xc == "rpbe":
		pp = "pbe"
	elif xc == "pw91":
		pp = "pw91"
	elif xc == "lda":
		pp = "lda"
	else:
		print("xc error")

	calc = EMT()

#	calc = Vasp(	prec=prec, xc=xc, pp=pp, ispin=2,
#					ismear=1, sigma=0.2,
#					isif=3,
#					ibrion=2, nsw=100, potim=potim, ediffg=ediffg,
#	    			kpts=kpts
#				)

	bulk.set_calculator(calc)
	bulk.get_potential_energy()

	a = bulk.cell[0,0] # optimized lattice constant

	return a

