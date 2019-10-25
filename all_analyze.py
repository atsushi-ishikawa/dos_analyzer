import os, json
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
		os.system('python analyze_dos.py ' + str(numpeaks) + ' ' + name)

json_to_csv(jsonfile, "tmp.csv")
#
# now edit csv file
#
df = pd.read_csv("tmp.csv")

# drop unnecessary features
del df["Unnamed: 0"];	del df["calculator"]
del df["constraints"];	del df["ctime"];		del df["dipole"];	del df["energy"];	del df["forces"]
del df["mtime"];		del df["numbers"];		del df["pbc"];		del df["positions"]
del df["tags"];			del df["unique_id"];	del df["user"]
del df["cell.array"];	del df["cell.pbc"];		del df["cell.__ase_objtype__"];			del df["lattice"]

if "stress" in df.columns:
	del df["stress"]

# format adjustment
orbital = ["s","p","d"]
keys    = ["height","position","width"]

for i in orbital:
	for j in keys:
		key = i + "-" + j
		tmp = df[key].str.replace("[","").str.replace("]","").str.split(",", expand=True)
		del df[key]

		for peak in range(numpeaks):
			df[i + '-' + j + str(peak+1)] = tmp[peak]

# convert to float
for i in df.columns:
	if i!="system":
		df[i] = df[i].astype(float)

df.set_index("system")

# delete system tag
del df["system"]

df = df.fillna(0.0)

# dropping strange values
print("dropping strange values...before:%4d" % len(df))

for i in orbital:
	for j in ["-height","-width"]:
		for k in range(1,numpeaks+1):
			key = i + j + str(k)
			if k==1:
				df = df[df[key] >  0.0] # zero is not allowed in the first peak -- at least one peak!
			else:
				df = df[df[key] >= 0.0]

print("dropping strange values...after: %4d" % len(df))

i = 0
col = df.iloc[:,i]
ave = np.mean(col)
std = np.std(col)
outlier_max = ave + 2*std
outlier_min = ave - 2*std

df = df[(df.iloc[:,i] < outlier_max) & (df.iloc[:,i] > outlier_min)]

df.dropna(how='any', axis=0)
print("dropping outliers...after: %4d"% len(df))

df = df[df["E_ads"] <  0.0]
print("dropping positive adsorption energy...after: %4d" % len(df))

df.to_csv(csvfile)
os.system("rm tmp.csv")

