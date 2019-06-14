import os
import pandas as pd
import numpy as np
from tools import json_to_csv

jsonfile = "tmpout.json"
csvfile  = "tmpout.csv"

numpeaks = 1
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
		if numpeaks==2:
			df[i + '-' + j + '2'] = tmp[1]

del df["Unnamed: 0"];  del df["calculator"]; del df["cell"]
del df["constraints"]; del df["ctime"]; del df["dipole"];    del df["energy"]; del df["forces"]
del df["lattice"];     del df["tags"];  del df["unique_id"]; del df["user"];   del df["pbc"]
del df["positions"];   del df["stress"]
del df["mtime"];       del df["numbers"]

# dropping strange values
print("dropping strange values...before:%4d" % len(df))
df = df[df["s-height1"].astype(float) > 0]
df = df[df["p-height1"].astype(float) > 0]
df = df[df["d-height1"].astype(float) > 0]

df = df[df["s-width1"].astype(float) > 0]
df = df[df["p-width1"].astype(float) > 0]
df = df[df["d-width1"].astype(float) > 0]

if numpeaks==2:
	df = df[df["s-height2"].astype(float) > 0]
	df = df[df["p-height2"].astype(float) > 0]
	df = df[df["d-height2"].astype(float) > 0]

	df = df[df["s-width2"].astype(float) > 0]
	df = df[df["p-width2"].astype(float) > 0]
	df = df[df["d-width2"].astype(float) > 0]

df.set_index("system")
#del df["system"]

print("dropping strange values...after: %4d" % len(df))

print("dropping outliers for Eads...before:%4d"% len(df))

i = 0
col = df.iloc[:,i]
ave = np.mean(col)
std = np.std(col)
outlier_max = ave + 2*std
outlier_min = ave - 2*std

df = df[(df.iloc[:,i] < outlier_max) & (df.iloc[:,i] > outlier_min)]

df.dropna(how='any', axis=0)
print("dropping outliers...after: %4d"% len(df))

df.to_csv(csvfile)
os.system("rm tmp.csv")

