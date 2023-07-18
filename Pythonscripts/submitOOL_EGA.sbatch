#!/bin/bash
#SBATCH --job-name=OOL_EGA
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=2G
#SBATCH --output=OOL_EGA.out
#SBATCH --error=OOL_EGA.err
#SBATCH -p normal
#SBATCH -c 4

module load python/3.9.0

python3 OOL_EGA.py
