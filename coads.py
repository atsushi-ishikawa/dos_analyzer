from ase import Atoms, Atom
from ase.calculators.vasp import Vasp
from ase.constraints import FixAtoms
from ase.optimize import BFGS
from ase.build import fcc111,add_adsorbate

surf = fcc111("Pt", size=[3,3,5], a=4.0, vacuum=15.0)
c = FixAtoms(indices=[atom.index for atom in surf if atom.tag >= 3])
surf.set_constraint(c)

calc_surf = Vasp(prec="normal", xc="PBE", ispin=2, algo="VeryFast",
		 encut=400.0, ibrion=2, nsw=100, istart=0,
		 kpts=[3,3,1],
		 lreal=True, npar=12, nsim=12
	   )

surf.set_calculator(calc_surf)
e_surf = surf.get_potential_energy()
#
#
#
cell = [10.0, 10.0, 10.0]
mol  = Atoms("CO", positions=[(0,0,0),(0,0,1.2)], cell=cell)
calc_mol = Vasp(prec="normal", xc="PBE", ispin=1, algo="VeryFast",
		encut=400.0, ibrion=2, nsw=100, istart=0,
		kpts=[1,1,1],
		lreal=True, npar=12, nsim=12
	   )

mol.set_calculator(calc_mol)

e_mol = mol.get_potential_energy()

add_adsorbate(surf, mol, 1.5, "ontop")

e_tot = surf.get_potential_energy()

e_ads = e_tot - (e_surf + e_mol)

print "Adsorption energy:", e_ads

