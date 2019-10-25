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
adsorbate = "O"
#adsorbate = "CO"
#adsorbate = "CH3"

if adsorbate=="CO":
	ads_height = 1.6
	ads_geom  = [(0, 0, 0), (0, 0, 1.2)]
elif adsorbate=="CH3":
	ads_height = 2.2
	ads_geom  = [(1, 1, 0), (0.4, 1, 0), (1.6, 1, 0), (1, 1.6, 0)]
else:
	ads_height = 1.8
	ads_geom  = [(0, 0, 0)]

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
	xc     = "pbe"
	prec   = "normal"
	encut  =  350
	nelmin =  5
	potim  =  0.10
	nsw    =  20
	ediff  =  1.0e-5
	ediffg = -0.03
	kpts   =  [2,2,1]
	gamma  =  True
	isym   = -1
	ispin  =  1 #### NOTICE: "analyze.dos" is not yet adjusted to ispin=2
	ibrion =  2
	nfree  =  20
	ispin_adsorbate = 2
	#
	# single point 
	#
	xc_sp     = "pbe"
	encut_sp  = 400
	ismear_sp = -5
	sigma_sp  = 0.1
	kpts_sp  = [3,3,1]

	npar = 18 # 18 for ito
	nsim = 18
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
surf_json = "surf_data.json"
# ads_json  = "ads_data.json"

surf_json = os.path.join(cudir, surf_json)
# ads_json  = os.path.join(cudir, ads_json)

db_surf   = connect(surf_json)

#
# Set calculator here. Different levels for optimization and single point (DOS)
#
if "vasp" in calculator:
	calc_surf_opt = Vasp(prec=prec, xc=xc, pp=pp, ispin=ispin, algo="VeryFast", 
						 encut=encut, ismear=1, sigma=0.2, istart=0, nelmin=nelmin, 
						 isym=isym, ibrion=ibrion, nsw=nsw, potim=potim, ediff=ediff, ediffg=ediffg,
			 			 kpts=kpts, gamma=gamma, npar=npar, nsim=nsim, lreal=True, nfree=nfree)
	calc_surf_sp  = Vasp(prec=prec, xc=xc_sp, pp=pp, ispin=ispin, algo="VeryFast", 
					     encut=encut_sp, ismear=ismear_sp, sigma=0.2, istart=0, nelmin=nelmin, 
					     isym=isym, ibrion=-1, nsw=0, potim=0, ediff=ediff, ediffg=ediffg,
			 		     kpts=kpts_sp, gamma=gamma, npar=npar, nsim=nsim, lreal=True, nfree=nfree, lorbit=10 )
elif "emt" in calculator:
	calc_surf_opt = EMT()
	calc_surf_sp  = calc_surf_opt
#
# ------------------------ bulk ---------------------------
#
if alloy:
	bulk = make_bulk(element1, element2=element2, comp1=comp1, repeat=repeat_bulk)
else:
	bulk = make_bulk(element, repeat=repeat_bulk)

## lattice optimization
lattice, a0 = lattice_info_guess(bulk)
#
print("old cell: ", bulk.cell)
optimize_lattice_constant(bulk, lattice=lattice, a0=a0, xc="PBEsol", 
						  encut=encut, ediff=ediff, ediffg=ediff*0.1, npar=npar, nsim=nsim)
print("new cell: ", bulk.cell)
quit()
#
# ------------------------ surface ------------------------
#
# load lattice constant form bulk calculaiton database
#
# surface construction
#
cell = bulk.get_cell()
bulk.set_cell(optimized_cell)
quit()
#print "cell,before",cell
#cell = cell // repeat_bulk
#print "cell,after", cell

surf = surface(bulk, face, nlayer, vacuum=vacuum)
surf.translate([0, 0, -vacuum])

surf = sort_atoms_by_z(surf)

calc_formation_energy = True
if alloy and calc_formation_energy:
	nat = 0
	e_bulk = 0.0
	for ielem in [element1, element2]:
		bulk = make_bulk(ielem, repeat=2)
		#a1 = get_optimized_lattice_constant(bulk, lattice=lattice, a0=a0, xc="PBEsol", 
		#									 encut=encut, ediff=ediff, ediffg=ediff*0.1, npar=npar, nsim=nsim)
		#
		# bulk energy for formation energy --- same with alloy surface calculator
		#

		# optimization
		bulk.set_calculator(calc_surf_opt)
		bulk.get_potential_energy()

		# single point
		bulk.set_calculator(calc_surf_sp)
		ene = bulk.get_potential_energy()

		ene = ene // len(bulk)
		nat = surf.get_chemical_symbols().count(ielem)
		e_bulk += ene * nat
		print("ene:{}".format(ene))
		print("nat:{}".format(nat))

	print("bulk_energy:", e_bulk)
#
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
# optimization
surf.set_calculator(calc_surf_opt)
surf.get_potential_energy()

# single point
surf.set_calculator(calc_surf_sp)
e_surf = surf.get_potential_energy()

print("surface number:",len(surf))
print("e_surf:{}".format(e_surf))
print("formation energy:", e_surf-e_bulk)

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
	efermi = calc_surf_sp.read_fermi()
	print("fermi energy:", efermi)
#
# ------------------------ adsorbate ------------------------
#
cell = [10.0, 10.0, 10.0]
mol  = Atoms(adsorbate, positions=ads_geom, cell=cell)

if "vasp" in calculator:
	calc_mol  = Vasp(prec=prec, xc=xc, pp=pp, ispin=ispin_adsorbate, algo="VeryFast",
					 encut=encut_sp, ismear=0, sigma=0.05, istart=0, nelmin=nelmin, 
					 isym=isym, ibrion=2, nsw=nsw, potim=potim, ediff=ediff, ediffg=ediffg,
					 kpts=[1,1,1], gamma=gamma, npar=npar, nsim=nsim, lreal=True, nfree=nfree)
elif "emt" in calculator:
	calc_mol = EMT()

mol.set_calculator(calc_mol)
e_mol = mol.get_potential_energy()
#
# ------------------- surface + adsorbate -------------------
#
if position_str == "atop":
	position = (0, 0)
	offset = (0.5, 0.5)
elif position_str == "hcp":
	position = (0,0)
	offset = (0.1667, 0.1667)
elif position_str == "fcc":
	position = (0,0)
	offset = (0.3333, 0.3333)

# shift H of CH3 to upper
if adsorbate=="CH3":
	for i,j in enumerate(mol.get_chemical_symbols()):
		if j=="H":
			mol[i].position[2] += 0.5

add_adsorbate(surf, mol, ads_height, position=position, offset=offset)
#
# optimization
#
surf.get_potential_energy()
#
# single point
#
if "vasp" in calculator:
	calc_surf_sp = Vasp(prec=prec, xc=xc_sp, pp=pp, ispin=ispin, algo="VeryFast", 
					    encut=encut_sp, ismear=ismear_sp, sigma=0.2, istart=0, nelmin=nelmin, 
					    isym=isym, ibrion=-1, nsw=0, potim=0, ediff=ediff, ediffg=ediffg,
			 		    kpts=kpts_sp, gamma=gamma, npar=npar, nsim=nsim, lreal=True, nfree=nfree)
elif "emt" in calculator:
	calc_surf_sp = EMT()

surf.set_calculator(calc_surf_sp)
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
			  #data={ adsorbate + "-" + position_str: e_ads} )
			  data={ "E_ads" : e_ads} )
#
# remove working directory
#
shutil.rmtree(workdir)

