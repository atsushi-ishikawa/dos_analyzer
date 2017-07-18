from vasptools import *
from ase.calculators.vasp import VaspDos
import matplotlib.pylab as plt

dos = VaspDos(doscar="DOSCAR_Pd111")
ene = dos.energy
newdos = smear_dos(dos)

plt.plot(ene,newdos)
plt.show()
