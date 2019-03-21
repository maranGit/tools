#!/bin/sh
#
# -------------------------------------------------------------------- #
#                                                                      #
#        bash script to submit 1 job on Habanero using slurm           #
#                                                                      #
#                   last modified: 3/17/2019 RM                        #
#                                                                      #
# -------------------------------------------------------------------- #
#
# The account name for the job.
#SBATCH --account=sun
#
# The job name
#SBATCH --job-name=polycrystal_3d_post
#
# The number of cpu cores to use
#SBATCH --cpus-per-task=1
#
# wall clock time, hour:minutes:seconds
#SBATCH --time=01:00:00
#
# The memory the job will use per cpu core
#SBATCH --mem-per-cpu=64gb
#
# User notification
#SBATCH --mail-type=NONE
#SBATCH --mail-user=rm3681@columbia.edu
#
# Error and log displayed on terminal
#SBATCH -o simulation.%j.log
#SBATCH -e simulation.%j.err
# 
# bin folder, store frequently used links
# 
PATH=$PATH:~/bin
pat2exii_3d <temp >post.out
exit
 
# End of script
