#!/bin/bash
#$ -cwd
#$ -j y
#$ -l h_data=16G,h_rt=8:00:00
#$ -e ./OutFiles
#$ -o ./OutFiles

. /u/local/Modules/default/init/modules.sh
module use $HOME/modulefiles
module load localpy
module load python/3.1


python3 variants_discoverer.py ../../data/chr1/hw2grad_M_1_chr_1.txt ../../data/chr1/reads_hw2grad_M_1_chr_1.txt

