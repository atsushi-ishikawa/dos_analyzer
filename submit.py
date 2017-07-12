import os,sys

dict = {
  "Rh" : {"lattice" : "fcc" , "a0" : 3.803}, 
  "Pd" : {"lattice" : "fcc" , "a0" : 3.891}, 
  "Ir" : {"lattice" : "fcc" , "a0" : 3.839}, 
  "Pt" : {"lattice" : "fcc" , "a0" : 3.924}
}

for element in dict:
	lattice = dict[element]["lattice"]
	a0      = dict[element]["a0"]
	print("element:%s, lattice:%s, a0=%f" % (element, lattice, a0))
	command = "qsub run.sh " + str(element) + " " + str(lattice) + " " + str(a0)
	os.system(command)

