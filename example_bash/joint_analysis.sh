#!/bin/bash

#SBATCH -N 1 
#SBATCH -p pe2,bigmem,dev
#SBATCH --mem=50G
#SBATCH -t 0-12:00 # Runtime in D-HH:MM
#SBATCH -J gruyere_joint # <-- name of job

#load required modules
conda activate gruyere
python ../src/gruyere_joint.py ../example_data/inputs.yaml 
