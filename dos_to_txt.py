from ase import Atoms, Atom
from ase.calculators.vasp import VaspDos
from vasptools import smear_dos
import sys
import numpy as np
import matplotlib.pylab as plt
import seaborn as sb

argvs = sys.argv
system  = argvs[1]
#system : "Pd_111"

if len(argvs) == 3:
	orbital = argvs[2]
	draw_pdos = True
else:
	draw_pdos = False

doscar = "DOSCAR_" + system
sigma = 10.0

#
# finding natom
#
f = open(doscar, "r")
line1 = f.readline()
natom = int( line1.split()[0] )
f.close()

dos = VaspDos(doscar=doscar)

ene  = dos.energy
tdos = dos.dos
tdos = smear_dos(ene, tdos, sigma=sigma)

outname = system + ".txt"
fout = open(outname,"w")

for i,x in enumerate(ene):
	fout.write("{0:<12.4f}{1:12.4e}\n".format(ene[i],tdos[i]))

fout.close()
