import glob, os
from ase import Atoms, Atom
from ase.io import read, write

#l = glob.glob("beef*/CONTCAR")

#for i, contcar in enumerate(l):
#	if "gas" in contcar:
#		continue
#	mol = read(contcar)
#	mol.set_constraint()
#	mol.wrap(eps=0.1)
#	write("CONTCAR", mol)
#	os.system('ase gui --rotations 280x,350y,359z CONTCAR --repeat 1,1,1 --output %.3d.png' % i)

for file in os.listdir():
	b,e = os.path.splitext(file)
	if e == ".xml":
		mol = read(b+e)
		mol.wrap(eps=0.1)
		write("CONTCAR", mol)
		#os.system("ase gui --rotations 290x,0y --repeat 3,2,1 --image-number -1 CONTCAR --output %s.png" % b)
		os.system("ase gui --rotations 0x,0y --repeat 1,1,1 --image-number -1 CONTCAR --output %s.png" % b)
