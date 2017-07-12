import os

files = ["WAVECAR","OSZICAR","XDATCAR","EIGENVAL","PCDAT","IBZKPT"]

for file in files:
	os.system("rm -f %s" % file)
