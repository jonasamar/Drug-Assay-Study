#!/bin/bash
#SBATCH --job-name=OOL_EGA_pen_other_models
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=2G
#SBATCH --output=OOL_EGA_pen_other_models.out
#SBATCH --error=OOL_EGA_pen_other_models.err
#SBATCH -p normal
#SBATCH -c 4

module load python/3.9.0

python3 OOL_EGA_pen_other_models.py