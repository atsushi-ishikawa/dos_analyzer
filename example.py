from dos_analyzer.vaspdosplus import VaspDosPlus
import os
import glob
from tinydb import TinyDB

db = TinyDB("sample.json")

doscardir = "doscars"
doscars = glob.glob(doscardir + "/" + "DOSCAR*")

descriptors = {}
for idoscar in doscars:
	system = idoscar.split("_")[1] + "_" + idoscar.split("_")[2]
	print(system)
	dos = VaspDosPlus(doscar=idoscar, system=system, numpeaks=1)
	dos.load_surface_data(json="surf_data.json")
	descriptor = dos.get_descriptors()
	db.insert(descriptor)
