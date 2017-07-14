import os,sys
#
# lattice constants takne from 
# http://periodictable.com/Properties/A/LatticeConstants.html
#
dict = {
  "Rh" : {"lattice" : "fcc" , "a0" : 3.803}, 
  "Pd" : {"lattice" : "fcc" , "a0" : 3.891}, 
  "Ir" : {"lattice" : "fcc" , "a0" : 3.839}, 
  "Pt" : {"lattice" : "fcc" , "a0" : 3.924},
  "Cu" : {"lattice" : "fcc" , "a0" : 3.615}, 
  "Ag" : {"lattice" : "fcc" , "a0" : 4.085}, 
  "Au" : {"lattice" : "fcc" , "a0" : 4.078} 
}

for element in dict:
	lattice = dict[element]["lattice"]
	a0      = dict[element]["a0"]
	command = "qsub run.sh " + str(element) + " " + str(lattice) + " " + str(a0)
	os.system(command)

