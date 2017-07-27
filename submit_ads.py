import os,sys

dict = {
  "Rh" : {"lattice" : "fcc" },
  "Pd" : {"lattice" : "fcc" }, 
  "Ir" : {"lattice" : "fcc" }, 
  "Pt" : {"lattice" : "fcc" },
  "Cu" : {"lattice" : "fcc" }, 
  "Ag" : {"lattice" : "fcc" }, 
  "Au" : {"lattice" : "fcc" } 
}

for element in dict:
	command = "qsub run_ads.sh " + str(element)
	os.system(command)

