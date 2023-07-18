#!/bin/bash
#SBATCH --job-name=corrEGA_noEGA_pen
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=2G
#SBATCH --output=corrEGA_noEGA_pen.out
#SBATCH --error=corrEGA_noEGA_pen.err
#SBATCH -p normal
#SBATCH -c 4

module load python/3.9.0

python3 corrEGA_noEGA_pen.py