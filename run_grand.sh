#!/bin/sh 
#------ pjsub option --------#
#PJM -L rscgrp=n22240a
#PJM -L node=1
#PJM --mpi proc=40
#PJM -L elapse=20:00:00 
#PJM -g n22240
#PJM -j
#------- Program execution -------#
module load intel

PRG=/lustre0/home/n22240/vasp/vasp.5.4.4/bin/vasp_std
OUT=stdout_$$

VASP_VDW=/lustre0/home/n22240/vasp/vasp.5.4.4/vdw_kernel.bindat

cp $VASP_VDW ./

# remove old results
rm stdout* stderr*

# INP variables are given by submit_ads.py: INP1 = element1, INP2 = adsorbate

#python adsorption.py ${INP1} ${INP2} ${INP3} 1> std_${INP1}_${INP2}_${INP3}_ads.out 2> err_${INP1}_${INP2}_${INP3}_ads.out # alloy
#python adsorption.py ${INP1} ${INP2} ${INP3} 1> std_${INP1}_${INP2}_${INP3}_ads.out 2> err_${INP1}_${INP2}_${INP3}_ads.out # alloy
#python adsorption.py --calculator="vasp" 1> stdout 2> stderr
#python adsorption.py --element1=${INP1} --adsorbate=${INP2} --calculator="vasp" 1> stdout 2> stderr
python adsorption.py --element1=${INP1} --adsorbate=${INP2} --calculator="vasp" 1> std_${INP1}_${INP2}.out 2> err_${INP1}_${INP2}.out

