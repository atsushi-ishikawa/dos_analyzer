#!/bin/bash
#PJM -L "rscunit=ito-a"
#PJM -L "rscgrp=ito-a-oc170117"
#PJM -L "vnode=1"
#PJM -L "vnode-core=36"
#PJM -L "elapse=10:00:00"
#PJM -j
#PJM -X

NUM_PROCS=36

module load intel/2017

PRG=${HOME}/vasp/vasp.5.4.4/bin/vasp_std
MOL=$(basename $PWD)
OUT=$MOL.out

export I_MPI_PERHOST=$NUM_CORES
export I_MPI_FABRICS=shm:ofa

export I_MPI_HYDRA_BOOTSTRAP=rsh
export I_MPI_HYDRA_BOOTSTRAP_EXEC=/bin/pjrsh
export I_MPI_HYDRA_HOST_FILE=${PJM_O_NODEINF}

echo $INP1 $INP2 $INP3

# mpiexec.hydra -n $NUM_PROCS $PRG > $OUT
# python adsorption.py $INP 1> std_${INP}_ads.out 2> err_${INP}_ads.out
python adsorption.py ${INP} 1>std_${INP}_ads.out 2>err_${INP}_ads.out
# python adsorption.py ${INP1} ${INP2} ${INP3} 1> std_${INP1}_${INP2}_${INP3}_ads.out 2> err_${INP1}_${INP2}_${INP3}_ads.out
