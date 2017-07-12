from ase import Atoms, Atom
from ase.calculators.vasp import Vasp
from ase.constraints import FixAtoms
from ase.optimize import BFGS
from ase.build import fcc111,add_adsorbate

from ase.db import connect

#
# Determine lattice constant.
# Take from bulk database.
#
bulk_json = "bulk_data.json"

db = connect(bulk_json)

element = "Rh"
xc = "pbesol"

atoms = db.get_atoms(element=element,xc=xc)
a = atoms.cell[0,0]

size = [2,2,5]
# --- INCAR keywords ---
prec = "normal"
encut = 400.0
potim = 0.1
ediff  = 1.0e-4
ediffg = -0.01
kpts=[3, 3, 1]
# ----------------------

xc = xc.lower()
if xc == "pbe" or xc == "pbesol" or xc == "rpbe":
	pp = "pbe"
elif xc == "pw91":
	pp = "pw91"
elif xc == "lda":
	pp = "lda"
else:
	print("xc error")

surf = fcc111(element, size=size, a=a, vacuum=15.0)
c = FixAtoms(indices=[atom.index for atom in surf if atom.tag >= 3])
surf.set_constraint(c)

calc_surf = Vasp(	prec=prec, xc=xc, pp=pp, ispin=2, algo="VeryFast",
			encut=encut, ismear=1, sigma=0.2, istart=0,
			ibrion=2, nsw=100, potim=potim, ediffg=ediffg,
		    	kpts=kpts, npar=12, nsim=12, lreal=True
       		    )

surf.set_calculator(calc_surf)
e_surf = surf.get_potential_energy()
#
#
#
cell = [10.0, 10.0, 10.0]
mol  = Atoms("CO", positions=[(0,0,0),(0,0,1.2)], cell=cell)

calc_mol = Vasp(prec=prec, xc=xc, pp=pp, ispin=1, algo="VeryFast",
		encut=encut, ismear=0, istart=0,
		ibrion=2, nsw=100, potim=potim, ediffg=ediffg,
		kpts=[1,1,1], npar=12, nsim=12, lreal=True
       		)

mol.set_calculator(calc_mol)

e_mol = mol.get_potential_energy()

add_adsorbate(surf, mol, 1.6, "ontop")

e_tot = surf.get_potential_energy()

e_ads = e_tot - (e_surf + e_mol)

print "Adsorption energy:", e_ads

