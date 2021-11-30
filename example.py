from dos_analyzer.vaspdosplus import VaspDosPlus
import os
import glob
import argparse
from tinydb import TinyDB

parser = argparse.ArgumentParser()
parser.add_argument("--jsonfile", default="sample.json", help="name of the json file containing surface information")
parser.add_argument("--numpeaks", default=1, type=int, help="number of peaks")
args = parser.parse_args()

jsonfile = args.jsonfile
numpeaks = args.numpeaks

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
	dos.numpeaks = numpeaks

	descriptor = dos.get_descriptors()

	data.update(descriptor)
	db.insert(data)
