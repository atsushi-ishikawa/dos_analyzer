import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import pyplot as plt
from matplotlib.collections import PolyCollection
from ase.calculators.vasp import VaspDos
import glob
from vasptools import get_efermi_from_doscar

def drawPropagation(xmin, xmax, z, doslist):
	sy = z.size
	#
	# temporary read DOSCAR for getting size
	#
	dos = VaspDos(doscar=doslist[0])
	T   = dos.energy
	idx = np.where( (xmin < T) & (T < xmax) ) # limit xmin-xmax
	T   = T[idx]

	sx = T.size
	T = np.tile(T, (sy, 1)) # extend T
	z = np.tile(z, (sx, 1)).T

	for i,idoscar in enumerate(doslist):
		dos    = VaspDos(doscar=idoscar)
		efermi = get_efermi_from_doscar(idoscar)
		T[i]  += efermi

	U = dos.dos[idx]
	U = np.tile(U, (sy, 1))

	for i,idoscar in enumerate(doslist):
		dos    = VaspDos(doscar=idoscar)
		U[i]   = dos.dos[idx]
		U[i,0] = 0.0 ; U[i,-1] = 0.0 # to make bottom line flat
		
	fig = plt.figure(figsize=(12,12))
	ax  = fig.add_subplot(1,1,1, projection="3d", proj_type="ortho")

	verts = []
	for i,idoscar in enumerate(doslist):
		efermi = get_efermi_from_doscar(idoscar)
		verts.append(list(zip(T[i,:]-efermi, U[i,:]/np.max(U[i,:]))))

	poly = PolyCollection(verts, facecolors=(1, 1, 1, 1.0), edgecolors=(0, 0, 0, 1), linewidth=1.4) # RGBA

	ax.add_collection3d(poly, zs=z[:,0], zdir="y")
	ax.set_xlim3d(xmin, xmax)
	ax.set_ylim3d(np.min(z)-2.0, np.max(z)-2.0)
	ax.set_zlim3d(0, 1.1)

	labels = []
	for idoscar in range(len(doslist)):
		labels.append(doslist[idoscar].split("_")[1])

	ax.set_yticks(np.linspace(np.min(z), np.max(z), len(doslist)))
	ax.set_yticks(z[:,0])
	ax.set_yticklabels(labels, va='center', ha='left', fontsize=14, fontname="Arial")

	ax.grid(False)

	ax.set_xticks(np.arange(xmin, xmax+1, 2))
	ax.set_xticklabels(np.arange(xmin, xmax+1, 2), fontsize=14, fontname="Arial")
	ax.set_zticks([])

	ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
	ax.w_zaxis.set_pane_color((0.0, 0.0, 0.0, 1.0))

	ax.xaxis.pane.set_edgecolor('black')
	ax.yaxis.pane.set_edgecolor('black')
	ax.xaxis.pane.fill = False
	ax.yaxis.pane.fill = False
	ax.zaxis.pane.fill = False

	ax.get_proj = lambda: np.dot(Axes3D.get_proj(ax), np.diag([1.0, 5.0, 0.1, 1]))
	ax.view_init(elev=15, azim=270)

# select num doscar files, based on adsorption energy

from ase.db import connect
json = "surf_data.json"
db   = connect(json)
num  = 20

doslist = glob.glob("DOSCAR_*")
E_ads_list = {}
for name in doslist:
	system = name.split("_")[1] + "_" + name.split("_")[2]
	try:
		id = db.get(system=system).id
		E_ads = db.get(id=id).data["E_ads"]
		E_ads_list[system] = E_ads
	except:
		pass

E_ads_list = sorted(E_ads_list.items(), key=lambda x:x[1], reverse=True)

doslist = []
for i in range(num):
	doscar = "DOSCAR" + "_" + E_ads_list[i][0]
	print(E_ads_list[i])
	doslist.append(doscar)

drawPropagation(-10, 2, np.linspace(-1, 1, len(doslist)), doslist)

plt.subplots_adjust(left=0.125, right=0.9, top=0.9, bottom=0.5)
plt.savefig("tmp.eps")
#plt.show()
