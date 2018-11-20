from ase import Atoms, Atom
from ase.calculators.vasp import VaspDos
from vasptools import smear_dos
import sys
import numpy as np
import matplotlib.pylab as plt
import seaborn as sb

argvs = sys.argv
system  = argvs[1]
#system : "Pd111"

if len(argvs) == 3:
	orbital = argvs[2]
	draw_pdos = True
else:
	draw_pdos = False

doscar = "DOSCAR_" + system
sigma = 20.0

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

if draw_pdos:
	pdos = np.zeros(len(tdos))
	for i in range(0,natom):
		pdos = pdos + dos.site_dos(i, orbital)
	pdos = smear_dos(ene, pdos, sigma=sigma)

sb.set(context='notebook', style='darkgrid', palette='deep', font='sans-serif', font_scale=1, color_codes=False, rc=None)

#plt.plot(ene, tdos,"r-",linewidth=2)
plt.plot(ene, tdos)
filename = "DOS_" + system + ".png"
plt.ylabel("Density of state (-)")
plt.xlabel("Energy (eV)")
plt.savefig(filename)
#plt.show()
