from ase import Atoms, Atom
from ase.calculators.vasp import Vasp
from ase.constraints import FixAtoms
from ase.optimize import BFGS
from ase.build import bulk, surface, add_adsorbate
from ase.db import connect

from vasptools import *

import os,sys
import shutil
import numpy as np

from ase.visualize import view # debugging
#
# Determine lattice constant.
#
# basic conditions
#
argvs     = sys.argv
element1  = argvs[1]

face      = (1,1,1) ; face_str = ",".join( map(str,face) ).replace(",","")
adsorbate = "CO"
ads_geom  = [(0, 0, 0), (0, 0, 1.2)]
position  = (0,0) # ontop
#
# computational
#
xc     = "pw91"
vacuum = 10.0
nlayer = 3
nrelax = 2
repeat_bulk = 2
#
# INCAR keywords
#
prec   = "low"
encut  =  350.0
potim  =  0.1
nsw    =  0
ediff  =  1.0e-2
ediffg = -0.05
kpts   = [3, 3, 1]
#
# directry things
#
cudir   = os.getcwd()
# workdir = os.path.join(cudir, element, face_str, adsorbate)
# workdir = os.path.join(cudir, element1 + "_" + face_str + "_" + adsorbate)
workdir = os.path.join(cudir, "tmpdir")
# shutil.rmtree(workdir)
os.makedirs(workdir)
os.chdir(workdir)
#
# database
#
bulk_json = "bulk_data.json"
surf_json = "surf_data.json"
ads_json  = "ads_data.json"

bulk_json = os.path.join(cudir, bulk_json)
surf_json = os.path.join(cudir, surf_json)
ads_json  = os.path.join(cudir, ads_json)

db_bulk   = connect(bulk_json)
db_surf   = connect(surf_json)
db_ads    = connect(ads_json)
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
# ------------------------ bulk ---------------------------
#
bulk = make_bulk(element1, repeat=repeat_bulk)
lattice, a0 = lattice_info_guess(bulk)
a = get_optimized_lattice_constant(bulk, lattice=lattice, a0=a0)
#
# ------------------------ surface ------------------------
#
# load lattice constant form bulk calculaiton database
#
# surface construction
#
cell = [a, a, a]
bulk.set_cell(cell)
surf = surface(bulk, face, nlayer, vacuum=vacuum)
surf = sort_atoms_by_z(surf)
#
# setting tags for relax/freeze
#
natoms   = len(surf.get_atomic_numbers())
one_surf = natoms / nlayer / repeat_bulk
tag = np.ones(natoms, int)
for i in range(natoms-1, natoms-nrelax*one_surf-1, -1):
	tag[i] = 0

surf.set_tags(tag)
#
# making surface
#
c = FixAtoms(indices=[atom.index for atom in surf if atom.tag == 1])
surf.set_constraint(c)
#
# calulate
#
calc_surf = Vasp(	prec=prec, xc=xc, pp=pp, ispin=1, algo="VeryFast",
			encut=encut, ismear=1, sigma=0.2, istart=0,
			ibrion=2, nsw=nsw, potim=potim, ediffg=ediffg,
		    	kpts=kpts, npar=12, nsim=12, lreal=True, lorbit=10
       		    )
surf.set_calculator(calc_surf)
e_surf = surf.get_potential_energy()
#
# copy DOSCAR
#
dosfile  = "DOSCAR_" + element1 + face_str
dosfile  = os.path.join(cudir, dosfile)
os.system("cp DOSCAR %s" % dosfile)
efermi = calc_surf.read_fermi()

print "fermi energy:",efermi

# db_surf.write(surf, element=element1, lattice=lattice, face=face_str)
#
# ------------------------ adsorbate ------------------------
#
cell = [10.0, 10.0, 10.0]
mol  = Atoms(adsorbate, positions=ads_geom, cell=cell)

calc_mol = Vasp(prec=prec, xc=xc, pp=pp, ispin=1, algo="VeryFast",
		encut=encut, ismear=0, istart=0,
		ibrion=2, nsw=nsw, potim=potim, ediffg=ediffg,
		kpts=[1,1,1], npar=12, nsim=12, lreal=True, lorbit=10
       		)

mol.set_calculator(calc_mol)
e_mol = mol.get_potential_energy()
#
# ------------------- surface + adsorbate -------------------
#
add_adsorbate(surf, mol, 1.8, position=position)
#
e_tot = surf.get_potential_energy()
#
e_ads = e_tot - (e_surf + e_mol)
#
print "Adsorption energy:", e_ads

