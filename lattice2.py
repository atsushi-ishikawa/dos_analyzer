def get_optimized_lattice_constant(element, laticce="fcc", xc="PBEsol"):
	""" functiona to return optimized bulk constant
	"""
	from ase import Atoms
	from ase.build import bulk
	from ase.calculators.vasp import Vasp
	from ase.optimize import BFGS
	import os,sys
	import shutil
	#
	# compuational condition for Vasp
	#
	prec   = "normal"
	potim  = 0.1
	ediff  = 1.0e-6
	ediffg = -0.001
	kpts = [10, 10, 10]
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
	if element in dict:
		print "I know this."
		lattice = dict[element]["lattice"]
		a0 = dict[element]["a0"]
	else:
		print "I don't know this."

	"""
	bulk = bulk(element, lattice, a=a0, cubic=True)

	xc = xc.lower()
	if xc == "pbe" or xc == "pbesol" or xc == "rpbe":
		pp = "pbe"
	elif xc == "pw91":
		pp = "pw91"
	elif xc == "lda":
		pp = "lda"
	else:
		print("xc error")

	calc = Vasp(	prec=prec, xc=xc, pp=pp, ispin=2,
			ismear=1, sigma=0.2,
			isif=3,
			ibrion=2, nsw=100, potim=potim, ediffg=ediffg,
	    	kpts=kpts
       		    )
	
	bulk.set_calculator(calc)
	bulk.get_potential_energy()

	a = bulk.cell[0,0] # optimized lattice constant

	return a
 	
	"""
