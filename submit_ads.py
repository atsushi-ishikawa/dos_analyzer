import os
import sys

# no need to do separate calculation for pure metals when making alloy-including database,

add_pure_metals = True  # set true for the first time
alloy = False
env = "grand"  # "ito" or "grand"

dict = {
  "0":  {"element": "Rh", "lattice": "fcc"},
  "1":  {"element": "Pd", "lattice": "fcc"},
  "2":  {"element": "Ir", "lattice": "fcc"},
  "3":  {"element": "Pt", "lattice": "fcc"},
  "4":  {"element": "Cu", "lattice": "fcc"},
  "5":  {"element": "Ag", "lattice": "fcc"},
  "6":  {"element": "Au", "lattice": "fcc"},
  "7":  {"element": "Ru", "lattice": "fcc"},
  "8":  {"element": "Os", "lattice": "fcc"},
  "9":  {"element": "Cd", "lattice": "fcc"},
  "10": {"element": "Hg", "lattice": "fcc"}
  }

if alloy:
    for i in range(0, len(dict)):
        elm1 = dict[str(i)]["element"]
        for j in range(i+1, len(dict)):
            elm2 = dict[str(j)]["element"]

            # --- change by 50% --> gives 21 systems w.o. pure metals
            #for cmp1 in range(50,100,50): # start, end, diff

            # --- change by 25% --> gives 63 systems w.o. pure metals
            #for cmp1 in range(25,100,25): # start, end, diff

            # --- change by 10%
            #first_time = False
            #for cmp1 in range(10,50,10):  # start, end, diff # part1
            #for cmp1 in range(50,100,10): # start, end, diff # part2

            # --- change by 5%
            first_time = True
            # for cmp1 in range(25,50,5):  # start, end, diff # part2
            # for cmp1 in range(50,75,5):  # start, end, diff # part3
            # for cmp1 in range(75,100,5): # start, end, diff # part4

            for cmp1 in range(5, 25, 5):   # start, end, diff # part1
                if env == "ito":
                    cmd = "pjsub -x \"INP1={0}\" -x \"INP2={1}\" -x \"INP3={2}\" run_ito.sh".format(elm1, elm2, cmp1)
                elif env == "grand":
                    cmd = "pjsub -x \"INP1={0}\" -x \"INP2={1}\" -x \"INP3={2}\" run_grand.sh".format(elm1, elm2, cmp1)

                os.system(command)

# do pure metal if the loop is first_time
adsorbates = ["CO", "CH3", "NO"]

if add_pure_metals:
    for i in range(0, len(dict)):
        elm = dict[str(i)]["element"]
        for adsorbate in adsorbates:
            if env == "ito":
                command = "pjsub run_ito.sh   -x \"INP1={0}\" \"INP2={1}\"".format(elm, adsorbate)
            elif env == "grand":
                command = "pjsub run_grand.sh -x INP1={0} -x INP2={1}".format(elm, adsorbate)
            os.system(command)
