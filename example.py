from dos_analyzer.vaspdosplus import VaspDosPlus
import os
import glob
from tinydb import TinyDB

jsonfile = "sample.json"
os.system("rm {}".format(jsonfile))
db = TinyDB(jsonfile)

doscardir = "doscars"
doscars = glob.glob(doscardir + "/" + "DOSCAR*")

for idoscar in doscars:
	system = idoscar.split("_")[1] + "_" + idoscar.split("_")[2]
	data = {"system": system}

	dos = VaspDosPlus(doscar=idoscar, system=system, numpeaks=2)
	dos.load_surface_data(json="surf_data.json")

	descriptor = dos.get_descriptors()

	data.update(descriptor)
	db.insert(data)
