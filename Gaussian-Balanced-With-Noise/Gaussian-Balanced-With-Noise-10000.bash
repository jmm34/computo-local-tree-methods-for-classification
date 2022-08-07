#!/bin/bash
#SBATCH -t 6-0
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 40
#SBATCH -J gaussian-balanced-with-noise
R CMD BATCH --no-restaure --no-save Example-Gaussian-Balanced-With-Noise-10000.R
