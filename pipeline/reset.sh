#!/bin/bash

#SBATCH -p main
#SBATCH -J reset
#SBATCH -t 1:00:00
#SBATCH -c 8
#SBATCH --mem=24G

rm -rf cage
rm -rf txrevise
rm -rf annots
