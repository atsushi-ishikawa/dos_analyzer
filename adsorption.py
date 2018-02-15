#
# usage: python [this_script] Pt (single metal)
#        python [this_script] Pt Pd 50 (alloy)
#
from ase import Atoms, Atom
from ase.calculators.vasp import Vasp
from ase.constraints import FixAtoms
from ase.build import bulk, surface, add_adsorbate
from ase.db import connect
from ase.visualize import view # debugging

from vasptools import *

import os, sys, shutil
import numpy as np

#
# --- Determine lattice constant ---
#

#
# basic conditions
#
argvs     = sys.argv
element1  = argvs[1]

if len(argvs) == 4:
	#
	# in case of alloys
	#
	alloy = True
	element2 = argvs[2]
	comp1    = int(argvs[3])
	comp2    = 100 - comp1
	element  = element1 + "{:.1f}".format(comp1/100.0) + element2 + "{:.1f}".format(comp2/100.0)
else:
	alloy = False
	element  = element1

face      = (1,1,1) ; face_str = ",".join( map(str,face) ).replace(",","")
adsorbate = "CO"
ads_geom  = [(0, 0, 0), (0, 0, 1.2)] # CO
# adsorbate = "CH3"
# ads_geom  = [(0, 0, 0), (-0.6, 0, 1.1), (0.6, 0, 1.1), (0, 0.6, 1.1)]
position  = (0,0)  ; position_str = "ontop" 

#
# computational
#
xc     = "pbe"
vacuum = 8.0
nlayer = 3
nrelax = 2
repeat_bulk = 2
#
# INCAR keywords
#
prec   = "low"
encut  =  350.0
potim  =  0.1
nsw    =  100
ediff  =  1.0e-4
ediffg = -0.03
kpts   = [2, 2, 1]
ispin  = 1 #### NOTICE: "analyze.dos" is not yet adjusted to ispin=2
#
# directry things
#
cudir   = os.getcwd()
workdir = os.path.join(cudir, element + "_" + face_str + "_" + adsorbate)
os.makedirs(workdir)
os.chdir(workdir)
#
# database to save data
#
# surf_json = "surf_data.json"
surf_json = "surf_data_alloy.json"
# ads_json  = "ads_data.json"

surf_json = os.path.join(cudir, surf_json)
# ads_json  = os.path.join(cudir, ads_json)

db_surf   = connect(surf_json)
# db_ads    = connect(ads_json)
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
if alloy:
	bulk = make_bulk(element1, element2=element2, comp1=comp1, repeat=repeat_bulk)
else:
	bulk = make_bulk(element, repeat=repeat_bulk)

bulk_copy = bulk

lattice, a0 = lattice_info_guess(bulk)
a = get_optimized_lattice_constant(bulk, lattice=lattice, a0=a0)

# a = 4.0 * repeat_bulk
# a = 7.80311
# a = a/repeat_bulk

# make bulk again
#if alloy:
#	bulk = make_bulk(element1, element2=element2, comp1=comp1, a0=a, repeat=repeat_bulk)
#else:
#	bulk = make_bulk(element, a0=a, repeat=repeat_bulk)

#
# ------------------------ surface ------------------------
#
# load lattice constant form bulk calculaiton database
#
# surface construction
#

cell = bulk.get_cell()
#print "cell,before",cell
#cell = cell/repeat_bulk
#print "cell,after", cell

bulk.set_cell(cell)
surf = surface(bulk, face, nlayer, vacuum=vacuum)
surf.translate([0, 0, -vacuum])

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
# ==================================================
#                calulation starts here
# ==================================================
#
calc_surf = Vasp(prec=prec, xc=xc, pp=pp, ispin=ispin, algo="VeryFast",
		 encut=encut, ismear=1, sigma=0.2, istart=0,
		 ibrion=2, nsw=nsw, potim=potim, ediffg=ediffg,
		 kpts=kpts, npar=12, nsim=12, lreal=True, lorbit=10 )
surf.set_calculator(calc_surf)
e_surf = surf.get_potential_energy()
#
# copy DOSCAR
#
dosfile  = "DOSCAR_" + element + "_" + face_str
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

calc_mol = Vasp(prec=prec, xc=xc, pp=pp, ispin=ispin, algo="VeryFast",
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

system = element + "_" + face_str
db_surf.write(surf, system=system, lattice=lattice,
			  data={ adsorbate + "-" + position_str: e_ads} )
#
# remove working directory
#

# shutil.rmtree(workdir)
