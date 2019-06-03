import os
import pandas as pd
from tools import json_to_csv

jsonfile = "tmpout.json"
csvfile  = "tmpout.csv"

os.system('rm %s' % jsonfile)

name_list = os.listdir("./")
for name in name_list:
	if "DOSCAR" in name:
		#name = name.replace("DOSCAR_","")
		name = name + " s p d"
		os.system('python analize_dos.py' + ' ' + name)

json_to_csv(jsonfile, "tmp.csv")

# now edit csv file
df = pd.read_csv("tmp.csv")

list1 = ["s","p","d"]
list2 = ["height","position","width"]

for i in list1:
	for j in list2:
		key = i + "-dos " + j
		tmp = df[key].str.replace("[","").str.replace("]","").str.split(",", expand=True)
		del df[key]

		df[i + '-' + j + '1'] = tmp[0]
		#df[i + '-' + j + '2'] = tmp[1]

del df["Unnamed: 0"];  del df["calculator"]; del df["cell"]
del df["constraints"]; del df["ctime"]; del df["dipole"];    del df["energy"]; del df["forces"]
del df["lattice"];     del df["tags"];  del df["unique_id"]; del df["user"];   del df["pbc"]
del df["positions"];   del df["stress"]
del df["mtime"];       del df["numbers"]

df.to_csv(csvfile)
os.system("rm tmp.csv")

