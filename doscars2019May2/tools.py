def ABcoord(mol,A,B):
	from ase import Atoms
	import numpy as np

	symbols = np.array( mol.get_chemical_symbols() )
	A_idx   = np.where(symbols==A)[0]
	B_list  = np.where(symbols==B)[0]
	AB_dist = mol.get_distances(A_idx, B_list)

	R_AB = np.min(AB_dist)
	coordinatingB = B_list[np.argmin(AB_dist)]

	return R_AB,coordinatingB

def run_packmol(xyz_file, a, num, outfile):
	import os

	packmol = "/Users/ishi/packmol/packmol"
	filetype = "xyz"

	cell1 = [0.0, 0.0, 0.0, a, a, a]
	cell2 = " ".join(map(str, cell1))

	f = open("pack_tmp.inp", "w")
	text = [
		"tolerance 2.0"             + "\n",
		"output "     + outfile     + "\n",
		"filetype "   + filetype    + "\n",
		"structure "  + xyz_file    + "\n",
		"  number "   + str(num)    + "\n",
		"  inside box " + cell2     + "\n",
		"end structure"
		]
	f.writelines(text)
	f.close()

	run_string = packmol + " < pack_tmp.inp"

	os.system(run_string)

	# os.system("rm pack_tmp.inp")

def json_to_csv(jsonfile, csvfile):
	import json
	import pandas as pd
	from pandas.io.json import json_normalize
	f = open(jsonfile, "r")
	d = json.load(f)

	dd = []
	nrec = len(d)
	for i in range(1, nrec):
		if str(i) in d:
			tmp = d[str(i)]
			dd.append(json_normalize(tmp))

	ddd = pd.concat(dd)

	newcol = []
	for key in ddd.columns:
		key = key.replace("calculator_parameters.", "")
		key = key.replace("key_value_pairs.", "")
		key = key.replace("data.", "")
		newcol.append(key)

	ddd.columns = newcol

	# sort data by "num"
	if "num" in ddd.columns:
		ddd2 = ddd.set_index("num")
		ddd  = ddd2.sort_index()

	ddd.to_csv(csvfile)

def load_ase_json(jsonfile):
	import json
	import pandas as pd
	from pandas.io.json import json_normalize
	f = open(jsonfile, "r")
	d = json.load(f)

	dd = []
	nrec = len(d)
	for i in range(1, nrec):
		if str(i) in d:
			tmp = d[str(i)]
			dd.append(json_normalize(tmp))

	ddd = pd.concat(dd)

	newcol = []
	for key in ddd.columns:
		key = key.replace("calculator_parameters.", "")
		key = key.replace("key_value_pairs.", "")
		key = key.replace("data.", "")
		newcol.append(key)

	ddd.columns = newcol

	# sort data by "num"
	if "num" in ddd.columns:
 		ddd2 = ddd.set_index("num")
 		ddd  = ddd2.sort_index()

	return ddd

def delete_num_from_json(num,jsonfile):
	from ase.db import connect
	import sys

	db = connect(jsonfile)
	id = db.get(num=num).id
	db.delete([id])

