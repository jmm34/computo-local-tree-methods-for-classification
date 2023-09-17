#!/bin/bash
#SBATCH -t 6-0
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 40
#SBATCH -J population-genetics
R CMD BATCH --no-restaure --no-save Example-Population-Genetics.R 
