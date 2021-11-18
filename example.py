from dos_analyzer.vaspdosplus import VaspDosPlus
import os
import glob
from tinydb import TinyDB

db = TinyDB("sample.json")

doscardir = "doscars"
doscars = glob.glob(doscardir + "/" + "DOSCAR*")

descriptors = {}
for idoscar in doscars:
	#print(idoscar)
	dos = VaspDosPlus(doscar=idoscar)
	descriptor = dos.get_descriptors()
	db.insert(descriptor)
