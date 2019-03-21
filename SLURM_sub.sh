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
#SBATCH --job-name=polycrystal_3d
#
# The number of cpu cores to use
#SBATCH --cpus-per-task=9
#
# wall clock time, hour:minutes:seconds
#SBATCH --time=07:00:00
#
# The memory the job will use per cpu core
#SBATCH --mem-per-cpu=2gb
#
# User notification
#SBATCH --mail-type=ALL
#SBATCH --mail-user=rm3681@columbia.edu
#
# Error and log displayed on terminal
#SBATCH -o simulation.%j.log
#SBATCH -e simulation.%j.err
# 
# bin folder, store frequently used links
# 
PATH=$PATH:~/bin
#
#   Reset the modules used by code.
#
source /etc/profile.d/modules.sh
module load intel-parallel-studio/2019
export INTEL=/rigel/opt/parallel_studio_xe_2019/compilers_and_libraries_2019.0.117/linux
#
#   Make directory called ${WORKDIR} on the compute node
#
SERVER=$SLURM_SUBMIT_HOST
WORKDIR=/rigel/sun/users/${USER}/SLURM_$SLURM_JOB_ID
SCP=/usr/bin/scp
SSH=/usr/bin/ssh
mkdir -p $WORKDIR

PERMDIR=${SLURM_SUBMIT_DIR}

SERVPERMDIR=${SLURM_SUBMIT_HOST}:${PERMDIR}

###############################################################
#                                                             #
#    Transfer files from server to local disk.                #
#                                                             #
###############################################################

stagein()
{

cd ${SLURM_SUBMIT_DIR}
cp -u * ${WORKDIR}
cd ${WORKDIR}

}

############################################################
#                                                          #
#    Execute the run.  Do not run in the background.       #
#                                                          #
############################################################

runprogram()
{

### Define number of processors
NPROCS=$SLURM_CPUS_PER_TASK

### Run a parallel OpenMP executable.
export OMP_NUM_THREADS=$NPROCS
FFT_finite_3d.exe <JC_test.inp >${SLURM_SUBMIT_DIR}/woutput

}

###########################################################
#                                                         #
#   Copy necessary files back to permanent directory.     #
#                                                         #
###########################################################

stageout()
{

 cp -u * ${SLURM_SUBMIT_DIR}

}

##################################################
#                                                #
#   Staging in, running the job, and staging out #
#   were specified above as functions.  Now      #
#   call these functions to perform the actual   #
#   file transfers and program execution.        #
#                                                #
##################################################

stagein
runprogram
# stageout 

###############################################################
#                                                             #
#   The epilogue script automatically deletes the directory   #
#   created on the local disk (including all files contained  #
#   therein.                                                  #
#                                                             #
###############################################################

# cd ${WORKDIR}
# rm * 
# cd $PBS_O_WORKDIR

exit
 
# End of script
