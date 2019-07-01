#
# usage: python [this_script] Pt (single metal)
#        python [this_script] Pt Pd 50 (alloy)
#
from ase import Atoms, Atom
from ase.calculators.vasp import Vasp
from ase.calculators.emt import EMT
from ase.constraints import FixAtoms
from ase.build import bulk, surface, add_adsorbate, niggli_reduce
from ase.db import connect
from ase.visualize import view # debugging
from ase.io import read,write

from vasptools import *

import os, sys, shutil
import numpy as np

calculator = "vasp"; calculator = calculator.lower()
#
# --- Determine lattice constant ---
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
	element  = element1 + "{:.2f}".format(comp1/100.0) + element2 + "{:.2f}".format(comp2/100.0)
else:
	alloy = False
	element  = element1

face = (1,1,1) ; face_str = ",".join( map(str,face) ).replace(",","")

position_str = "atop" # atop, hcp, fcc
#adsorbate = "O"
adsorbate = "CO"
#adsorbate = "CH3"

#ads_geom  = [(0, 0, 0)]
ads_geom  = [(0, 0, 0), (0, 0, 1.2)]
#ads_geom  = [(0, 0, 0), (-0.6, 0, 1.1), (0.6, 0, 1.1), (0, 0.6, 1.1)]

ads_height = 1.6

vacuum = 10.0
nlayer = 2
nrelax = 2
repeat_bulk = 2
#
# computational
#
if "vasp" in calculator:
	#
	# INCAR keywords
	#
	xc     = "rpbe"
	prec   = "normal"
	encut  =  400
	nelmin =  5
	potim  =  0.10
	nsw    =  200
	ediff  =  1.0e-6
	ediffg = -0.03 # -0.03
	kpts   = [3,3,1]
	gamma  = True
	isym   = 0
	ispin  = 1 #### NOTICE: "analyze.dos" is not yet adjusted to ispin=2
	ibrion = 1
	nfree  = 20
	ispin_adsorbate = 1

	npar = 40 # 18 for ito
	nsim = 40
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
		pp = "pbe"

	## --- EMT --- -> nothing to set

#
# directry things
#
cudir   = os.getcwd()
workdir = os.path.join(cudir, element + "_" + face_str + "_" + adsorbate)
#workdir = os.path.join("/tmp/" + element + "_" + face_str + "_" + adsorbate) # whisky
os.makedirs(workdir)
os.chdir(workdir)
shutil.copy("../vdw_kernel.bindat" , ".")
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
# ------------------------ bulk ---------------------------
#
if alloy:
	bulk = make_bulk(element1, element2=element2, comp1=comp1, repeat=repeat_bulk)
else:
	bulk = make_bulk(element, repeat=repeat_bulk)

#bulk_copy = bulk

### pos 0
#a = 7.80
### pos 0 --> good

## lattice optimization
### pos 1
lattice, a0 = lattice_info_guess(bulk)
a = get_optimized_lattice_constant(bulk, lattice=lattice, a0=a0)
### pos 1 --> fail

# bulk_copy = make_bulk(element1, element2=element2, comp1=comp1, a0=a, repeat=repeat_bulk)
# del bulk
# bulk = bulk_copy

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
#cell = cell // repeat_bulk
#print "cell,after", cell

bulk.set_cell(cell)
surf = surface(bulk, face, nlayer, vacuum=vacuum)
surf.translate([0, 0, -vacuum])

surf = sort_atoms_by_z(surf)
# setting tags for relax/freeze
#
natoms   = len(surf.get_atomic_numbers())
one_surf = natoms // nlayer // repeat_bulk
tag = np.ones(natoms, int)
for i in range(natoms-1, natoms-nrelax*one_surf-1, -1):
	tag[i] = 0

surf.set_tags(tag)
#
# if x goes to negative ... invert to positive
#
if np.any( surf.cell[:,0] < 0 ):
	#
	# invert atoms
	#
 	for i in range(len(surf)):
 		xpos = surf[i].position[0]
 		surf[i].position[0] = -xpos
	#
 	# invert cell
	#
 	xpos = surf.cell[1,0]
 	surf.cell[1,0] = -xpos

surf.wrap()
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
if "vasp" in calculator:
	calc_surf = Vasp(prec=prec, xc=xc, pp=pp, ispin=ispin, algo="VeryFast", 
					 encut=encut, ismear=1, sigma=0.2, istart=0, nelmin=nelmin, 
					 isym=isym,ibrion=ibrion, nfree=nfree, nsw=nsw, potim=potim, ediffg=ediffg,
			 		 kpts=kpts, gamma=gamma, npar=npar, nsim=nsim, lreal=True, 
					 lorbit=10 )
elif "emt" in calculator:
	calc_surf = EMT()

surf.set_calculator(calc_surf)
e_surf = surf.get_potential_energy()

if "vasp" in calculator:
	#
	# copy DOSCAR
	#
	dosfile  = "DOSCAR_" + element + "_" + face_str
	dosfile  = os.path.join(cudir, dosfile)
	os.system("cp DOSCAR %s" % dosfile)
	#
	# printout Efermi
	#
	efermi = calc_surf.read_fermi()
	print("fermi energy:",efermi)

# db_surf.write(surf, element=element1, lattice=lattice, face=face_str)
#
# ------------------------ adsorbate ------------------------
#
cell = [10.0, 10.0, 10.0]
mol  = Atoms(adsorbate, positions=ads_geom, cell=cell)

if "vasp" in calculator:
	calc_mol  = Vasp(prec=prec, xc=xc, pp=pp, ispin=ispin_adsorbate, algo="VeryFast",
					 encut=encut, ismear=0, sigma=0.05, istart=0, nelmin=nelmin, 
					 isym=isym, ibrion=ibrion, nfree=nfree, nsw=nsw, potim=potim, ediffg=ediffg,
					 kpts=[1,1,1], gamma=gamma, npar=npar, nsim=nsim, lreal=True, lorbit=10 )
elif "emt" in calculator:
	calc_mol = EMT()

mol.set_calculator(calc_mol)
e_mol = mol.get_potential_energy()
#
# ------------------- surface + adsorbate -------------------
#
if position_str == "atop":
	# position = (a*2.0/11.0, a*1.0/11.0) # when nlayer = 2
	position = (0,0) # when nlayer = 1
	offset = (0.5, 0.5)
elif position_str == "hcp":
	position = (0,0) # when nlayer = 1
	offset = (0.1667, 0.1667)
elif position_str == "fcc":
	position = (0,0) # when nlayer = 1
	offset = (0.3333, 0.3333)

add_adsorbate(surf, mol, ads_height, position=position, offset=offset)
#
e_tot = surf.get_potential_energy()
e_ads = e_tot - (e_surf + e_mol)
#
print("Adsorption energy:", e_ads)
#
# copy vasprun.xml
#
if "vasp" in calculator:
	xmlfile  = "vasprun_" + element + "_" + face_str + ".xml"
	xmlfile  = os.path.join(cudir, xmlfile)
	os.system("cp vasprun.xml %s" % xmlfile)
#
# write to surf
#
system = element + "_" + face_str
db_surf.write(surf, system=system, lattice=lattice,
			  data={ adsorbate + "-" + position_str: e_ads} )
#
# remove working directory
#
shutil.rmtree(workdir)

