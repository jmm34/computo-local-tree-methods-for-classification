#!/bin/bash
#SBATCH -t 6-0
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 40
#SBATCH -J gaussian-balanced-no-noise
R CMD BATCH --no-restaure --no-save Example-Gaussian-Balanced-No-Noise-10000.R
