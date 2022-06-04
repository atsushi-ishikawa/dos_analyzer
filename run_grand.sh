#!/bin/sh 
#------ pjsub option --------#
#PJM -L rscgrp=n22238a
#PJM -L node=1
#PJM --mpi proc=20
#PJM -L elapse=40:00:00 
#PJM -g n22238
#PJM -j
#------- Program execution -------#

module load intel

PRG=/lustre0/home/n22240/vasp/vasp.5.4.4/bin/vasp_std
OUT=stdout_$$

VASP_VDW=/lustre0/home/n22240/vasp/vasp.5.4.4/vdw_kernel.bindat

cp $VASP_VDW ./

# remove old results
rm stdout* stderr*

#mpiexec.hydra -n ${PJM_MPI_PROC} ${PRG} >& ${OUT}
#python adsorption.py ${INP} 1> std_${INP}_ads.out 2> err_${INP}_ads.out # pure metal
python adsorption.py ${INP1} ${INP2} ${INP3} 1> std_${INP1}_${INP2}_${INP3}_ads.out 2> err_${INP1}_${INP2}_${INP3}_ads.out # alloy

