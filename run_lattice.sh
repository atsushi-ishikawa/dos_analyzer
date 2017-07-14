#!/bin/bash

#$ -S /bin/bash
#$ -cwd
#$ -q all.q
#$ -pe openmpi12 12

element=$1
lattice=$2
a0=$3

python lattice.py $element $lattice $a0 1> std_$element.out 2> err_$element.out

