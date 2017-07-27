import os,sys

#dict = {
#  "Rh" : {"lattice" : "fcc" },
#  "Pd" : {"lattice" : "fcc" }, 
#  "Ir" : {"lattice" : "fcc" }, 
#  "Pt" : {"lattice" : "fcc" },
#  "Cu" : {"lattice" : "fcc" }, 
#  "Ag" : {"lattice" : "fcc" }, 
#  "Au" : {"lattice" : "fcc" } 
#}

element1 = "Rh"
element2 = "Pd"
for comp1 in range(25,100,25):
	command = "qsub run_alloy.sh " + str(element1) + " " + str(element2) + " " + str(comp1)
	print command
#for element in dict:
#	command = "qsub run_ads.sh " + str(element)
#	os.system(command)

