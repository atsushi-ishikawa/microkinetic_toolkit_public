#!/bin/bash

#$ -S /bin/bash
#$ -cwd
#$ -q all.q
#$ -pe openmpi12 12


python reaction_energy.py 1> stdout.txt 2> stderr.txt

