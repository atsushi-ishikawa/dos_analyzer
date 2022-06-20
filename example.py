from dos_analyzer.vaspdosplus import VaspDosPlus
import os
import glob
import argparse
from tinydb import TinyDB

parser = argparse.ArgumentParser()
parser.add_argument("--surf_json", default="surf_xy.json", help="name of the json file containing surface information")
parser.add_argument("--numpeaks",  default=1, type=int, help="number of peaks in surface dos")
args = parser.parse_args()

surf_json = args.surf_json
numpeaks = args.numpeaks

os.system("rm {}".format(surf_json))
db = TinyDB(surf_json)

doscardir = "."
doscars = glob.glob(doscardir + "/" + "DOSCAR*")

surface_doscars = []
for idos in doscars:
    if idos.count("_") == 2:
        surface_doscars.append(idos)

# adsorbate
adsorbates = ["CO", "CH3"]

data = {}
for adsorbate in adsorbates:
    dos = VaspDosPlus(doscar="DOSCAR" + "_" + adsorbate)
    dos.ads_json = adsorbate + "_x.json"
    dos.system = adsorbate
    dos.numpeaks = 1  # only one peak in adsorbate
    descriptor = dos.get_descriptors(adsorbate=True)
    data.update(descriptor)
    db.insert(data)

for idoscar in surface_doscars:
    surface = idoscar.split("_")[1] + "_" + idoscar.split("_")[2]
    data = {}

    dos = VaspDosPlus(doscar=idoscar)
    dos.surf_jsonfile = "surf_data.json"
    dos.system = surface + "_" + adsorbate
    dos.numpeaks = numpeaks
    descriptor = dos.get_descriptors()

    data.update(descriptor)
    db.insert(data)
