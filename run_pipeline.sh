#!/bin/bash
#SBATCH --job-name=muscle_hdr
#SBATCH --output=muscle_hdr_%j.out
#SBATCH --error=muscle_hdr_%j.err
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=4:00:00
#SBATCH --mem=16G

module load python/3.9
source /lustre/vishal.bharti/muscle_HDR_scRNA/muscle_env/bin/activate

cd /lustre/vishal.bharti/muscle_HDR_scRNA
snakemake --snakefile workflow/Snakefile --cores 4
