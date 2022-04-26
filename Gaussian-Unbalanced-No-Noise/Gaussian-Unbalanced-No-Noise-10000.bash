#!/bin/bash
#SBATCH -t 6-0
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 40
#SBATCH -J gaussian-unbalanced-no-noise
R CMD BATCH --no-restaure --no-save Example-Gaussian-Unbalanced-No-Noise-10000.R
