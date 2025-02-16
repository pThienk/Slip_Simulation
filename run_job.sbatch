#!/bin/bash

# Controlling options for SLURM
# !REMEMBER TO CHANGE TO YOUR OWN NETID!
##########################################################################################
#
#SBATCH --account=pt27                   		# Request owner account name
#SBATCH --time=04:00:00                  		# Job run time (hh:mm:ss)
#SBATCH --nodes=1                        		# Number of nodes
#SBATCH --ntasks-per-node=32             		# Number of task (cores/ppn) per node
#SBATCH --job-name=slip_sim_%A           		# Name of batch job
#SBATCH --partition=secondary            		# Partition (queue)
#SBATCH --output=scratch/sim_output_%A_%a_%j    # Name of batch job output file
#SBATCH --error=scratch/sim_error_%A_%a_%j      # Name of batch job error file
#SBATCH --mail-user=pt27@illinois.edu    		# Send email notifications
#SBATCH --mail-type=BEGIN,END            		# Type of email notifications to send
#
##########################################################################################

# Load the GCC toolchain (makes the compilers and library accessible)
module load gcc/13.3.0

# Load the Intel compilers module (Make the Intel compilers
# and libraries accessible in the user environment)
#module load intel/18.0


# Define the number of threads used (for OpenMP code only, defaults to 1)
#export OMP_NUM_THREADS=16

# Run OpenMP or Serial code
./slip_sim_cluster # Arguments go here

# Summary of batch job info saved to a file (Optional).
qstat -f $SLURM_JOB_ID > $SLURM_JOB_NAME-$SLURM_JOB_ID.info-summary
