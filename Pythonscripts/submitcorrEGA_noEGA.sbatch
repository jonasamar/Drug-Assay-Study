#!/bin/bash
#SBATCH --job-name=corrEGA_noEGA
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=2G
#SBATCH --output=corrEGA_noEGA.out
#SBATCH --error=corrEGA_noEGA.err
#SBATCH -p normal
#SBATCH -c 4

module load python/3.9.0

python3 corrEGA_noEGA.py
