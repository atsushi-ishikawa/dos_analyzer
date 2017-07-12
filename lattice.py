from ase import Atoms
from ase.build import bulk
from ase.calculators.vasp import Vasp
from ase.optimize import BFGS
from ase.constraints import StrainFilter
from ase.db import connect
import re
import sys

jsonfile = "bulk_data.json"

argvs = sys.argv

element = argvs[1]
lattice = argvs[2]
a0      = argvs[3]

prec = "normal"
potim = 0.1
ediff  = 1.0e-6
ediffg = -0.001
kpts=[10, 10, 10]

bulk = bulk(element, lattice, a=a0, cubic=True)
# bulk = bulk.repeat((2,2,2))

for xc in ["PBEsol", "PBE", "PW91", "RPBE"]:

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

	trajectory = element + "_" + "lattice.traj"

	print("optimized lattice constant : %f" % bulk.cell[0,0])

	#
	# what to write
	#
	db = connect(jsonfile)
	db.write(bulk,
		data = {"element" : element, "lattice" : lattice} )

