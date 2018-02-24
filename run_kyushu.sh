#!/bin/bash
#PJM -L "rscunit=ito-a"
#PJM -L "rscgrp=ito-a-oc170117"
#PJM -L "vnode=1"
#PJM -L "vnode-core=12"
#PJM -L "elapse=10:00"
#PJM -j
#PJM -X

NUM_PROCS=12

module load intel/2017

PRG=${HOME}/vasp/vasp.5.4.4/bin/vasp_std
MOL=`basename $PWD`
OUT=$MOL.out

export I_MPI_PERHOST=$NUM_CORES
export I_MPI_FABRICS=shm:ofa

export I_MPI_HYDRA_BOOTSTRAP=rsh
export I_MPI_HYDRA_BOOTSTRAP_EXEC=/bin/pjrsh
export I_MPI_HYDRA_HOST_FILE=${PJM_O_NODEINF}

# mpiexec.hydra -n $NUM_PROCS $PRG > $OUT
python adsorption.py ads.txt

