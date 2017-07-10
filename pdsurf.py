from ase import Atoms, Atom
from ase.calculators.vasp import Vasp
from ase.constraints import FixAtoms
from ase.optimize import BFGS
from ase.build import fcc111

surf = fcc111("Pd", size=[3,3,5], a=4.0, vacuum=15.0)
c = FixAtoms(indices=[atom.index for atom in surf if atom.tag >= 3])
surf.set_constraint(c)

calc = Vasp(	prec="normal", xc="PBE", ispin=2, algo="VeryFast",
		encut=400.0, ibrion=2, nsw=100,
		kpts=[3,3,1],
		lreal=True, npar=12, nsim=12
	   )

surf.set_calculator(calc)
surf.get_potential_energy()

# dyn=BFGS(surf,trajectory="tmp.traj")
# dyn.run(fmax=0.05)

