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
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--element1", default="Pt")
parser.add_argument("--element2", default=None)
parser.add_argument("--comp1", default=None)
parser.add_argument("--adsorbate", default="CH3")
parser.add_argument("--surf_json", default="surf_data.json")
parser.add_argument("--calculator", default="emt")

args = parser.parse_args()
element1   = args.element1
element2   = args.element2
comp1      = args.comp1
adsorbate  = args.adsorbate
surf_json  = args.surf_json
calculator = args.calculator.lower()

optimize_lattice = True
calc_formation_energy = True  # alculate formation energy of BULK ALLOY from its component metals
do_cohp = False
lobster = "/lustre0/home/n22240/lobster/lobster-3.2.0/lobster-3.2.0"

if element2 is not None:
    alloy   = True
    comp2   = 100 - comp1
    element = element1 + "{:.2f}".format(comp1/100.0) + element2 + "{:.2f}".format(comp2/100.0)
else:
    alloy   = False
    element = element1

face = (1, 1, 1)
face_str = ",".join(map(str, face)).replace(",", "")

position_str = "atop"  # atop, hcp, fcc

if adsorbate in ["CO", "NO", "CH"]:
    ads_height = 1.6
    ads_geom  = [(0, 0, 0), (0, 0, 1.2)]
elif adsorbate == "CH3":
    ads_height = 2.2
    ads_geom  = [(1, 1, 0), (0.3, 0.6, 0), (1.7, 0.6, 0), (1, 1.8, 0)]
else:
    # atom
    ads_height = 1.8
    ads_geom  = [(0, 0, 0)]

vacuum = 10.0
nlayer = 2  # 3
nrelax = 2
repeat_bulk = 2
#
# computational
#
if "vasp" in calculator:
    # INCAR keywords
    xc     =  "pbe"
    ismear =  1
    sigma  =  0.1
    prec   = "normal"
    encut  =  450
    nelmin =  5
    nelm   =  40
    potim  =  0.10
    nsw    =  5  # 200
    ediff  =  1.0e-5
    ediffg = -0.05
    kpts   =  [1, 1, 1]
    gamma  =  True
    isym   =  -1
    ispin  =  1  # NOTICE: "analyze.dos" is not yet adjusted to ispin=2
    ibrion =  2
    nfree  =  20
    ispin_adsorbate = 1
    lreal  =  True

    # single point
    xc_sp     = xc
    encut_sp  = encut
    ismear_sp = 1  # -5
    sigma_sp  = 0.1
    gamma_sp  = True
    kpts_sp   = [1, 1, 1]  # [4, 4, 1]

    npar = 10  # 18 for ito
    nsim = 10

    # xc set
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
#workdir = os.path.join("/tmp/" + element + "_" + face_str + "_" + adsorbate)  # whisky
workdir = os.path.join(cudir, element + "_" + face_str + "_" + adsorbate)  # other

clean = True
if os.path.isdir(workdir) and clean:
    shutil.rmtree(workdir)

os.makedirs(workdir)
os.chdir(workdir)
shutil.copy(os.environ["HOME"] + "/vasp/vasp.5.4.4/vdw_kernel.bindat", ".")
#
# database to save data
#
surf_json = os.path.join(cudir, surf_json)
db_surf   = connect(surf_json)
#
# Set calculator here. Different levels for optimization and single point (DOS)
#
if "vasp" in calculator:
    calc_opt = Vasp(prec=prec, xc=xc, pp=pp, ispin=ispin, algo="VeryFast",
                    encut=encut, ismear=ismear, sigma=sigma, istart=0, nelm=nelm, nelmin=nelmin,
                    isym=isym, ibrion=ibrion, nsw=nsw, potim=potim, ediff=ediff, ediffg=ediffg,
                    kpts=kpts, gamma=gamma, npar=npar, nsim=nsim, lreal=lreal, nfree=nfree)
    calc_sp  = Vasp(prec=prec, xc=xc_sp, pp=pp, ispin=ispin, algo="VeryFast",
                    encut=encut_sp, ismear=ismear_sp, sigma=sigma_sp, istart=0, nelm=nelm, nelmin=nelmin,
                    isym=isym, ibrion=-1, nsw=0, potim=0, ediff=ediff, ediffg=ediffg,
                    kpts=kpts_sp, gamma=gamma_sp, npar=npar, nsim=nsim, lreal=lreal, lorbit=10)
    kpts_bulk_sp = [max(kpts_sp) for i in range(3)]
    calc_bulk_sp = copy.deepcopy(calc_sp)
elif "emt" in calculator:
    calc_opt = EMT()
    calc_sp  = EMT()
    calc_bulk_sp = EMT()
#
# ------------------------ bulk ---------------------------
#
if alloy:
    bulk = make_bulk(element1, element2=element2, comp1=comp1, repeat=repeat_bulk)
else:
    bulk = make_bulk(element, repeat=repeat_bulk)

## lattice optimization
lattice, a0 = lattice_info_guess(bulk)
if "vasp" in calculator and optimize_lattice:
    print("calculating bulk --- lattice optimization", flush=True)
    bulk_opt = optimize_lattice_constant(bulk, lattice=lattice, a0=a0, xc=xc, encut=encut, isif=6,
                                         ediff=ediff, ediffg=ediff*0.1, nsw=nsw, npar=npar, nsim=nsim)
else:
    bulk_opt = bulk

bulk_opt.set_calculator(calc_bulk_sp)
e_bulk = bulk_opt.get_potential_energy()

e_form = 0.0
if alloy and calc_formation_energy:
    nat = 0
    e_bulk_component = 0.0
    for ielem in [element1, element2]:
        tmpbulk = make_bulk(ielem, repeat=2)
        bulk_opt = optimize_lattice_constant(tmpbulk, lattice=lattice, a0=a0, xc=xc, encut=encut,
                                             ediff=ediff, ediffg=ediff*0.1, nsw=nsw, npar=npar, nsim=nsim)

        # bulk energy for formation energy --- same with alloy surface calculator

        # single point
        bulk_optset_calculator(calc_bulk_sp)
        ene = bulk_opt.get_potential_energy()
        ene /= len(bulk_opt)
        nat = bulk_opt.get_chemical_symbols().count(ielem)
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
surf = surface(bulk_opt, face, nlayer, vacuum=vacuum, periodic=True)
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

    # invert atoms
    for i in range(len(surf)):
        xpos = surf[i].position[0]
        surf[i].position[0] = -xpos

    # invert cell
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
print("calculating surface --- optimization", flush=True)
surf.set_calculator(calc_opt)
surf.get_potential_energy()

# single point
print("calculating surface --- single point", flush=True)
if do_cohp:
    nbands = calc_opt.read_number_of_electrons()//2 + len(surf)*4  # needed to increase nbands for COHP by lobster
    calc_sp.int_params["nbands"] = nbands  # replace the calculator

surf.set_calculator(calc_sp)
e_slab = surf.get_potential_energy()

# surface energy
a, b, c, alpha, beta, gamma = surf.cell.cellpar()
surf_area = a*b*math.sin(math.radians(gamma))
e_surf = (e_slab - (len(surf)/len(bulk))*e_bulk) / (2.0*surf_area)

if "vasp" in calculator:
    # copy DOSCAR of surface to currentdir
    dosfile = "DOSCAR_" + element + "_" + face_str
    dosfile = os.path.join(cudir, dosfile)
    os.system("cp DOSCAR %s" % dosfile)

    # printout Efermi
    efermi = calc_sp.read_fermi()

    # do lobster and copy COHP file
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
mol  = Atoms(adsorbate, positions=ads_geom, cell=cell, pbc=True)

if "vasp" in calculator:
    calc_mol_opt  = Vasp(prec=prec, xc=xc, pp=pp, ispin=ispin_adsorbate, algo="VeryFast",
                         encut=encut_sp, ismear=0, sigma=0.05, istart=0, nelm=nelm, nelmin=nelmin,
                         isym=isym, ibrion=2, nsw=nsw, potim=potim, ediff=ediff, ediffg=ediffg,
                         kpts=[1, 1, 1], gamma=gamma, npar=npar, nsim=nsim, lreal=True, nfree=nfree)
    calc_mol_sp   = Vasp(prec=prec, xc=xc, pp=pp, ispin=ispin_adsorbate, algo="VeryFast",
                         encut=encut_sp, ismear=0, sigma=0.05, istart=0, nelm=nelm, nelmin=nelmin,
                         isym=isym, ibrion=-1, nsw=0, potim=potim, ediff=ediff, ediffg=ediffg,
                         kpts=[1, 1, 1], gamma=gamma, npar=npar, nsim=nsim, lreal=True, nfree=nfree, lorbit=10)
elif "emt" in calculator:
    calc_mol_opt = EMT()
    calc_mol_sp  = EMT()

print("calculating adsorbate --- {:s}".format(adsorbate), flush=True)

# opt
mol.set_calculator(calc_mol_opt)
mol.get_potential_energy()

# single point
mol.set_calculator(calc_mol_sp)
mol.get_potential_energy()

e_mol = mol.get_potential_energy()

# copy DOSCAR of adsorbate
dosfile = "DOSCAR_" + adsorbate
dosfile = os.path.join(cudir, dosfile)
if not os.path.exists(dosfile):
    os.system("cp DOSCAR %s" % dosfile)
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
print("calculating surface + adsorbate --- optimization", flush=True)
surf.set_calculator(calc_opt)
surf.get_potential_energy()
#
# single point
#
print("calculating surface + adsorbate --- single point", flush=True)
surf.set_calculator(calc_sp)
e_tot = surf.get_potential_energy()
e_ads = e_tot - (e_slab + e_mol)
#
# copy vasprun.xml
#
if "vasp" in calculator:
    xmlfile  = "vasprun_" + element + "_" + face_str + "_" + adsorbate + ".xml"
    xmlfile  = os.path.join(cudir, xmlfile)
    os.system("cp vasprun.xml %s" % xmlfile)
#
# write to surf
#
system = element + "_" + face_str + "_" + adsorbate
db_surf.write(surf, system=system, lattice=lattice, data={"E_ads": e_ads, "E_form": e_form, "E_surf": e_surf})
                                                          
print("done: E_ads = {:6.3f}".format(e_ads), flush=True)
