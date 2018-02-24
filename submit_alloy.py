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
		# for comp1 in range(0,150,50): # start, end, diff
		comp1 = 50
		print element1,element2,comp1
		# command = "qsub run_alloy.sh " + str(element1) + " " + str(element2) + " " + str(comp1) # gridengine
		command = "pjsub -x \"INP1={0}\" -x \"INP2={1}\" -x \"INP3={2}\" run_kyushu.sh".format(element1,element2,comp1) # PBS
		os.system(command)

