import os

name_list = os.listdir("./")
for name in name_list:
	if "DOSCAR" in name:
		name = name.replace("DOSCAR_","")
		name = name + " s p d"
		os.system('python analize_dos.py' + ' ' + name)
