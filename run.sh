#!/bin/bash

#$ -S /bin/bash
#$ -cwd
#$ -q all.q
#$ -pe openmpi12 12
#$ -l hostname=whisky??&!whisky04

INP="gri.txt"
LBL=$$

python reaction_energy.py $INP 1> stdout_$LBL.txt 2> stderr_$LBL.txt

