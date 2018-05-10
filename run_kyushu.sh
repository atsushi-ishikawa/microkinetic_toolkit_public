#!/bin/bash
#PJM -L "rscunit=ito-a"
#PJM -L "rscgrp=ito-a-oc170117"
#PJM -L "vnode=4"
#PJM -L "vnode-core=36"
#PJM -L "elapse=12:00:00"
#PJM -j
#PJM -X

NUM_PROCS=144

module load intel/2017

PRG=${HOME}/vasp/vasp.5.4.4/bin/vasp_std
MOL=`basename $PWD`
OUT=$MOL.out

echo "--------- $INP"

export I_MPI_PERHOST=$NUM_CORES
export I_MPI_FABRICS=shm:ofa

export I_MPI_HYDRA_BOOTSTRAP=rsh
export I_MPI_HYDRA_BOOTSTRAP_EXEC=/bin/pjrsh
export I_MPI_HYDRA_HOST_FILE=${PJM_O_NODEINF}

# mpiexec.hydra -n $NUM_PROCS $PRG > $OUT
# python adsorption.py ch4.txt
# python adsorption.py $INP

INP="test4.txt"
LBL=$$

####
vasp_script="/home/usr6/m70286a/ase/run_vasp.py"
echo 
echo "import os" > $vasp_script
echo "exitcode = os.system(\"mpirun -np ${NUM_PROCS} ${PRG}\")" >> $vasp_script
####

python reaction_energy.py $INP 1> stdout_$LBL.txt 2> stderr_$LBL.txt

