from ase import Atoms
from ase.build import bulk
from ase.calculators.vasp import Vasp
from ase.optimize import BFGS
from ase.constraints import StrainFilter

element = "Pt"
lattice = "fcc"
a0 = 3.4

potim = 0.1
ediffg = -0.01
kpts=[10,10,10]

bulk = bulk(element, lattice, a=a0, cubic=True)
# bulk = bulk.repeat((2,2,2))

calc = Vasp(	prec="normal", xc="PBE", 
		ismear=1, sigma=0.2,
		isif=3,
		ibrion=2, nsw=100, potim=potim, ediffg=ediffg,
	    	kpts=kpts
           )

bulk.set_calculator(calc)
bulk.get_potential_energy()

# sf = StrainFilter(bulk)

trajectory = element + "_" + "lattice.traj"

# opt = BFGS(bulk, trajectory=trajectory)
# opt.run(fmax=0.005)

print("optimized lattice constant : %f" % bulk.cell[0,0])

