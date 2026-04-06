#!/bin/bash

#SBATCH --partition=sapphire
#SBATCH --job-name="ell31"
#SBATCH --account="punim2890"
#SBATCH --ntasks=256
#SBATCH --mem-per-cpu=8000
#SBATCH --mail-user=m.kortge@unimelb.edu.au
#SBATCH --mail-type=ALL
#SBATCH --time=0-7:00:00

if [ -z "$SLURM_JOB_ID" ]; then
  echo "You need to submit your job to the queuing system with sbatch"
  exit 1
fi

cd /data/gpfs/projects/punim2890/class-group-tabulation || exit 1

module purge
module load gompi/2024a
module load GMP/6.3.0
module load Meson/1.4.0

meson setup --reconfigure build ./
meson compile -C build clgrp_ell ell_count

D_MAX=137438953472
FILES=512
ELL=31
FOLDER=../lmfdb
OUTPUT=./ell${ELL}.csv

srun -n "$SLURM_NTASKS" build/clgrp_ell "$D_MAX" "$FILES" "$ELL" "$FOLDER"
srun -n "$SLURM_NTASKS" build/ell_count "$D_MAX" "$FILES" "$ELL" "$FOLDER" >"$OUTPUT"

my-job-stats -a -n -s
