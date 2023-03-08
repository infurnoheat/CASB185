#!/bin/bash

qsub -cwd -N "cov_1_quant" -o $HOME/joblogs/csbio185/ -j y -l h_rt=24:00:00,h_data=8G -pe shared 8 -M $USER@mail -m bea -t 1-4:1 kallisto_quant.sh cov_1.txt 131 21.1
qsub -cwd -N "cov_2_quant" -o $HOME/joblogs/csbio185/ -j y -l h_rt=24:00:00,h_data=8G -pe shared 8 -M $USER@mail -m bea -t 1-4:1 kallisto_quant.sh cov_2.txt 128 22.2
qsub -cwd -N "moc_1_quant" -o $HOME/joblogs/csbio185/ -j y -l h_rt=24:00:00,h_data=8G -pe shared 8 -M $USER@mail -m bea -t 1-4:1 kallisto_quant.sh moc_1.txt 135 20.2
qsub -cwd -N "moc_2_quant" -o $HOME/joblogs/csbio185/ -j y -l h_rt=24:00:00,h_data=8G -pe shared 8 -M $USER@mail -m bea -t 1-4:1 kallisto_quant.sh moc_2.txt 134 20.6
