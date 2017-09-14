import os,sys

dict = {
  "0"  : {"element" : "Rh" , "lattice" : "fcc" },
  "1"  : {"element" : "Pd" , "lattice" : "fcc" }, 
  "2"  : {"element" : "Ir" , "lattice" : "fcc" }, 
  "3"  : {"element" : "Pt" , "lattice" : "fcc" }, 
  "4"  : {"element" : "Cu" , "lattice" : "fcc" }, 
  "5"  : {"element" : "Ag" , "lattice" : "fcc" }, 
  "6"  : {"element" : "Au" , "lattice" : "fcc" }
}

for i in range(0,len(dict)):
	element1 = dict[str(i)]["element"]
	for j in range(i+1,len(dict)):
		element2 = dict[str(j)]["element"]
		for comp1 in range(25,100,25):
			print element1,element2,comp1
			command = "qsub run_alloy.sh " + str(element1) + " " + str(element2) + " " + str(comp1)
			os.system(command)

