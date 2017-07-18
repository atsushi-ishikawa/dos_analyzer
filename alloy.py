def get_alloy(element1, element2, comp1=100,lattice="fcc",a0=4.0,repeat=2):
	from ase import Atoms
	from ase.build import bulk
	from ase.calculators.emt import EMT
	from ase.optimize import BFGS
	from ase.visualize import view
	import os,sys
	import shutil
	import numpy as np
	#
	# to fractional
	#
	comp2 = 100 - comp1
	#
	comp1 = comp1/100.0
	comp2 = comp2/100.0
	#
	# form bulk first
	#
	bulk = bulk(element1, lattice, a=a0, cubic=True)
	bulk = bulk.repeat( (repeat,repeat,repeat) )
	#
	#
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

