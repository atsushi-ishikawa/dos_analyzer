#!/bin/bash

#$ -S /bin/bash
#$ -cwd
#$ -q all.q
#$ -pe openmpi12 12

# python coads.py 1> std.out 2> err.out
python lattice.py Pt fcc 3.6 1> std.out 2> err.out

python clean.py
