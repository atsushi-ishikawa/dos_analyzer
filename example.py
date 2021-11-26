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
	data = {}

	dos = VaspDosPlus(doscar=idoscar)
	dos.surf_jsonfile = "surf_data.json"
	dos.system = system
	dos.numpeaks = 1

	descriptor = dos.get_descriptors()

	data.update(descriptor)
	db.insert(data)
