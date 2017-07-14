from ase import Atoms, Atom
from ase.calculators.vasp import Vasp
from ase.constraints import FixAtoms
from ase.optimize import BFGS
# from ase.build import fcc111,add_adsorbate
from ase.build import bulk, surface, add_adsorbate
from ase.db import connect

import os,sys
import shutil
import numpy as np

from ase.visualize import view # debugging
#
# Determine lattice constant.
# Take from bulk database.
#
bulk_json = "bulk_data.json"
db = connect(bulk_json)

#
# basic conditions
#
argvs = sys.argv
element   = argvs[1]
lattice   = "fcc"
face      = (1,1,1) ; face_str = ",".join( map(str,face) ).replace(",","")
adsorbate = "CO"
ads_geom  = [(0, 0, 0), (0, 0, 1.2)]
position  = "ontop"
#
# computational
#
xc = "pw91"
vacuum = 10.0
nlayer = 5
nrelax = 2
repeat_size = (3,3,1)
#
# load lattice constant form bulk calculaiton database
#
atoms = db.get_atoms(element=element, xc="pbesol")
a = atoms.cell[0,0]
#
# INCAR keywords
#
prec  = "low"
encut = 350.0
potim = 0.1
nsw   = 100
ediff  = 1.0e-4
ediffg = -0.05
kpts=[3, 3, 1]
#
# directry things
#
cudir   = os.getcwd()
# workdir = os.path.join(cudir, element, face_str, adsorbate)
workdir = os.path.join(cudir, element + "_" + face_str + "_" + adsorbate)
# shutil.rmtree(workdir)
os.makedirs(workdir)
os.chdir(workdir)
#
# xc set
#
xc = xc.lower()
if xc == "pbe" or xc == "pbesol" or xc == "rpbe":
	pp = "pbe"
elif xc == "pw91":
	pp = "pw91"
elif xc == "lda":
	pp = "lda"
else:
	print("xc error")
#
# surface construction
#
bulk = bulk(element, lattice, a=a, cubic=False)
surf = surface(bulk, face, nlayer, vacuum=vacuum)
#
# setting tags for relax/freeze
#
tag = np.ones(nlayer, int)
for i in range(nlayer-1, nlayer-nrelax-1, -1):
	tag[i] = 0
surf.set_tags(tag)

c = FixAtoms(indices=[atom.index for atom in surf if atom.tag == 1])
surf.set_constraint(c)
surf = surf.repeat(repeat_size)

calc_surf = Vasp(	prec=prec, xc=xc, pp=pp, ispin=1, algo="VeryFast",
			encut=encut, ismear=1, sigma=0.2, istart=0,
			ibrion=2, nsw=nsw, potim=potim, ediffg=ediffg,
		    	kpts=kpts, npar=12, nsim=12, lreal=True, lorbit=11
       		    )
surf.set_calculator(calc_surf)
e_surf = surf.get_potential_energy()
#
# copy DOSCAR
#
dosfile  = "DOSCAR_" + element + face_str
dosfile  = os.path.join(cudir, dosfile)
os.system("cp DOSCAR %s" % dosfile)
#
#
#
"""
cell = [10.0, 10.0, 10.0]
mol  = Atoms(adsorbate, positions=ads_geom, cell=cell)

calc_mol = Vasp(prec=prec, xc=xc, pp=pp, ispin=1, algo="VeryFast",
		encut=encut, ismear=0, istart=0,
		ibrion=2, nsw=nsw, potim=potim, ediffg=ediffg,
		kpts=[1,1,1], npar=12, nsim=12, lreal=True
       		)

mol.set_calculator(calc_mol)

e_mol = mol.get_potential_energy()

add_adsorbate(surf, mol, 1.6, position=position)

e_tot = surf.get_potential_energy()

e_ads = e_tot - (e_surf + e_mol)

print "Adsorption energy:", e_ads

"""
