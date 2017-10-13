#!/bin/bash

#$ -S /bin/bash
#$ -cwd
#$ -q all.q
#$ -pe openmpi12 12

element=$1

python adsorption.py $element 1> std_${element}_ads.out 2> err_${element}_ads.out

