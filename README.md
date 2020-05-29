## VASPでDOSを自動的に計算してmachine-learning (PCCP投稿中論文)

1. VASPの計算を行う:`adsorption.py`

### Plotting DOS

```bash {cmd="/bin/bash"}
scp whisky:/home/a_ishi/ase/ase_vasp/DOSCAR_Pd_111 ./DOSCAR_Pd111
```

```python {cmd="/Users/ishi/.pyenv/shims/python"}
from ase import Atoms, Atom
from ase.calculators.vasp import VaspDos
import sys
from vasptools import smear_dos
import numpy as np
import matplotlib.pylab as plt
import seaborn as sb

argvs = sys.argv
# system  = argvs[1]
system = "Pd111"

if len(argvs) == 3:
	orbital = argvs[2]
	draw_pdos = True
else:
	draw_pdos = False

doscar = "DOSCAR_" + system
sigma = 10.0

#
# finding natom
#
f = open(doscar, "r")
line1 = f.readline()
natom = int( line1.split()[0] )
f.close()

dos = VaspDos(doscar=doscar)

ene  = dos.energy
tdos = dos.dos
tdos = smear_dos(ene, tdos, sigma=sigma)

if draw_pdos:
	pdos = np.zeros(len(tdos))
	for i in range(0,natom):
		pdos = pdos + dos.site_dos(i, orbital)
	pdos = smear_dos(ene, pdos, sigma=sigma)

sb.set(context='notebook', style='darkgrid', palette='deep',
    font='sans-serif', font_scale=1, color_codes=False, rc=None)

#plt.plot(ene, tdos,"r-",linewidth=2)
plt.plot(ene, tdos)
filename = "DOS_" + system + ".png"
plt.ylabel("Density of state (-)")
plt.xlabel("Energy (eV)")
# plt.savefig(filename)
plt.show()
```

* to see original DOS
```bash
python plotdos.py Rh0.2Cu0.8_111 ; open DOS_Rh0.2Cu0.8_111.png
```

#### DOSのpeakを検出してjsonファイルに追記する
* 計算データの入ったjsonファイル("surf_data.json")が必要

##### 個別のDOSCARに行う場合
```bash
python analize_dos.py [system_name] [orbital(s/p/d)]
 example: python analize_dos.py Rh0.2Cu0.8_111 d
```

##### フォルダ内の全てのDOSCARに行う場合
```bash
python all_analize.py
```

#### After getting DOS
* convert to CSV
```python
from tools import json_to_csv
json_to_csv("asdf.json","asdf.csv")
```
* or analyze by pandas
```python
from tools import load_ase_json
df = load_ase_json("asdf.json") # dataframe
```

