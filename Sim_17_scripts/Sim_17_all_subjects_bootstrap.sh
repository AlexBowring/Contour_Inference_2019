#!/bin/bash
#$ -S /bin/bash
#$ -l h_vmem=8G
#$ -l h_rt=23:59:00
#$ -t 1
#$ -cwd
#$ -o $HOME/log
#$ -e $HOME/log

. /etc/profile

module add matlab

matlab -nodisplay -nojvm -r Sim_17_all_subjects_bootstrap

