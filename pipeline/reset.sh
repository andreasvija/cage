#!/bin/bash

#SBATCH -p main
#SBATCH -J reset
#SBATCH -t 1:00:00
#SBATCH -c 8
#SBATCH --mem=24G

rm -rf cage
rm -rf txrevise
rm -rf txrevise-25
rm -rf annots
rm -rf annots_5
rm -rf annots_10
rm -rf annots_15
rm -rf annots_20
rm -rf annots_25
rm -rf annots_50
rm -rf annots_100
rm -rf annots_150
rm -rf annots_200
