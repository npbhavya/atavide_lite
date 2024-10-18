#!/bin/bash

#SBATCH --job-name=atavide
#SBATCH --output=%x-%j.out.txt
#SBATCH --error=%x-%j.err.txt
#SBATCH --time=1-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=500G
#SBATCH --partition=high-capacity
#SBATCH --qos=hc-concurrent-jobs

#nextflow run atavide/workflow/install.nf 
nextflow run atavide/workflow/main.nf
