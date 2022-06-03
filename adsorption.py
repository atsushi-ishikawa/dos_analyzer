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
from ase.visualize import view  # debugging
from ase.io import read, write
from vasptools import *

import os
import sys
import shutil
import math
import copy
import numpy as np

calculator = "vasp"; calculator = calculator.lower()

# Whether to calculate formation energy of BULK ALLOY from its component metals.
calc_formation_energy = True
do_cohp = False
lobster = "/lustre0/home/n22240/lobster/lobster-3.2.0/lobster-3.2.0"
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

face = (1, 1, 1)
face_str = ",".join(map(str, face)).replace(",", "")

position_str = "atop"  # atop, hcp, fcc
#adsorbate = "O"
#adsorbate = "CO"
adsorbate = "CH3"

if adsorbate == "CO":
    ads_height = 1.6
    ads_geom  = [(0, 0, 0), (0, 0, 1.2)]
elif adsorbate == "CH3":
    ads_height = 2.2
    ads_geom  = [(1, 1, 0), (0.4, 1, 0), (1.6, 1, 0), (1, 1.6, 0)]
else:
    ads_height = 1.8
    ads_geom  = [(0, 0, 0)]

vacuum = 10.0
nlayer = 3
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
    ismear = 1
    sigma  = 0.1
    prec   = "normal"
    encut  =  400
    nelmin =  5
    nelm   =  100
    potim  =  0.10
    nsw    =  200
    ediff  =  1.0e-5
    ediffg = -0.05
    kpts   =  [2, 2, 1]
    gamma  =  True
    isym   =  -1
    ispin  =  1  # NOTICE: "analyze.dos" is not yet adjusted to ispin=2
    ibrion =  2
    nfree  =  20
    ispin_adsorbate = 2
    #
    # single point
    #
    xc_sp     = xc
    encut_sp  = encut
    ismear_sp = 1  # -5
    sigma_sp  = 0.1
    gamma_sp  = True
    kpts_sp   = [4, 4, 1]

    npar = 10  # 18 for ito
    nsim = 10
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
shutil.copy("../vdw_kernel.bindat", ".")
#
# database to save data
#
surf_json = "surf_data.json"
surf_json = os.path.join(cudir, surf_json)
db_surf   = connect(surf_json)
#
# Set calculator here. Different levels for optimization and single point (DOS)
#
if "vasp" in calculator:
    calc_opt = Vasp(prec=prec, xc=xc, pp=pp, ispin=ispin, algo="VeryFast",
                    encut=encut, ismear=ismear, sigma=sigma, istart=0, nelm=nelm, nelmin=nelmin,
                    isym=isym, ibrion=ibrion, nsw=nsw, potim=potim, ediff=ediff, ediffg=ediffg,
                    kpts=kpts, gamma=gamma, npar=npar, nsim=nsim, lreal=True, nfree=nfree)
    calc_sp  = Vasp(prec=prec, xc=xc_sp, pp=pp, ispin=ispin, algo="VeryFast",
                    encut=encut_sp, ismear=ismear_sp, sigma=sigma_sp, istart=0, nelm=nelm, nelmin=nelmin,
                    isym=isym, ibrion=-1, nsw=0, potim=0, ediff=ediff, ediffg=ediffg,
                    kpts=kpts_sp, gamma=gamma_sp, npar=npar, nsim=nsim, lreal=True, lorbit=10)
    kpts_bulk_sp = [max(kpts_sp) for i in range(3)]
    calc_bulk_sp = copy.deepcopy(calc_sp)
elif "emt" in calculator:
    calc_opt = EMT()
    calc_sp  = calc_opt
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
optimize_lattice_constant(bulk, lattice=lattice, a0=a0, xc=xc,
                          encut=encut, ediff=ediff, ediffg=ediff*0.1, npar=npar, nsim=nsim)
bulk.set_calculator(calc_bulk_sp)

e_bulk = bulk.get_potential_energy()

e_form = 0.0
if alloy and calc_formation_energy:
    nat = 0
    e_bulk_component = 0.0
    for ielem in [element1, element2]:
        tmpbulk = make_bulk(ielem, repeat=2)
        optimize_lattice_constant(tmpbulk, lattice=lattice, a0=a0, xc=xc,
                                  encut=encut, ediff=ediff, ediffg=ediff, npar=npar, nsim=nsim)
        #
        # bulk energy for formation energy --- same with alloy surface calculator
        #
        # single point
        tmpbulk.set_calculator(calc_bulk_sp)
        ene = tmpbulk.get_potential_energy()
        ene /= len(tmpbulk)
        #nat = surf.get_chemical_symbols().count(ielem)
        nat = bulk.get_chemical_symbols().count(ielem)
        e_bulk_component += ene * nat

    # formation energy of bulk alloy from its component
    e_form = e_bulk - e_bulk_component
#
# ------------------------ surface ------------------------
#
# load lattice constant form bulk calculaiton database
#
# surface construction
#
surf = surface(bulk, face, nlayer, vacuum=vacuum)
surf.translate([0, 0, -vacuum])

surf = sort_atoms_by_z(surf)

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
if np.any(surf.cell[:, 0] < 0):
    #
    # invert atoms
    #
    for i in range(len(surf)):
        xpos = surf[i].position[0]
        surf[i].position[0] = -xpos
    #
    # invert cell
    #
    xpos = surf.cell[1, 0]
    surf.cell[1, 0] = -xpos

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
surf.set_calculator(calc_opt)
surf.get_potential_energy()

# single point
nbands = calc_opt.read_number_of_electrons()//2 + len(surf)*4  # needed to increase nbands for COHP by lobster
calc_sp.int_params["nbands"] = nbands  # replace the calculator
surf.set_calculator(calc_sp)
e_slab = surf.get_potential_energy()

# surface energy
a, b, c, alpha, beta, gamma = surf.get_cell_lengths_and_angles()
surf_area = a*b*math.sin(math.radians(gamma))
e_surf = (e_slab - (len(surf)/len(bulk))*e_bulk) / (2.0*surf_area)

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
    efermi = calc_sp.read_fermi()
    #
    # do lobster and copy COHP file
    #
    if do_cohp:
        make_lobsterin()
        os.system("env OMP_NUM_THREADS=%d %s" % (npar, lobster))
        cohpfile  = "COHPCAR_" + element + "_" + face_str
        cohpfile  = os.path.join(cudir, cohpfile)
        os.system("cp COHPCAR.lobster %s" % cohpfile)
#
# ------------------------ adsorbate ------------------------
#
cell = [10.0, 10.0, 10.0]
mol  = ase.Atoms(adsorbate, positions=ads_geom, cell=cell)

if "vasp" in calculator:
    calc_mol  = Vasp(prec=prec, xc=xc, pp=pp, ispin=ispin_adsorbate, algo="VeryFast",
                     encut=encut_sp, ismear=0, sigma=0.05, istart=0, nelm=nelm, nelmin=nelmin,
                     isym=isym, ibrion=2, nsw=nsw, potim=potim, ediff=ediff, ediffg=ediffg,
                     kpts=[1, 1, 1], gamma=gamma, npar=npar, nsim=nsim, lreal=True, nfree=nfree)
elif "emt" in calculator:
    calc_mol = EMT()

mol.set_calculator(calc_mol)
e_mol = mol.get_potential_energy()
#
# ------------------- surface + adsorbate -------------------
#
if position_str == "atop":
    if nlayer == 2:
        position = (0, 0)
        offset = (0.5, 0.5)
    elif nlayer == 3:
        position = (0, 0)
        offset = (0.3333, 0.3333)
elif position_str == "hcp":
    position = (0, 0)
    offset = (0.1667, 0.1667)
elif position_str == "fcc":
    position = (0, 0)
    offset = (0.3333, 0.3333)

# shift H of CH3 to upper
if adsorbate == "CH3":
    for i, j in enumerate(mol.get_chemical_symbols()):
        if j == "H":
            mol[i].position[2] += 0.5

add_adsorbate(surf, mol, ads_height, position=position, offset=offset)
#
# optimization
#
surf.set_calculator(calc_opt)
surf.get_potential_energy()
#
# single point
#
surf.set_calculator(calc_sp)
e_tot = surf.get_potential_energy()
e_ads = e_tot - (e_slab + e_mol)
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
db_surf.write(surf, system=system, lattice=lattice, data={"E_ads": e_ads, "E_form": e_form, "E_surf": e_surf})
#
# remove working directory
#
shutil.rmtree(workdir)
