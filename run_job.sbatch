#!/bin/bash

# Controlling options for SLURM
# !REMEMBER TO CHANGE TO YOUR OWN NETID!
##########################################################################################
#
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

logfile="slip_batch_submission.log"
input_file="slip_sim_input.txt"  # Change this to your input file with numerical values
i=1
# Note that $(something) calls the 'something' command and replaces $(something) with the output.
echo "Starting batch submissions at $(date)" > $logfile

# Loop over each value in the input file
while IFS=',' read -r epsilon1 epsilon2; do
    # Skip empty lines or lines that start with # (commented lines)
    if [[ -z "$epsilon1" || -z "$epsilon2" || "$epsilon1" =~ ^# ]]; then
        continue
    fi

    # Export the two EPSILON values as environment variables
    export EPSILON1=$epsilon1
    export EPSILON2=$epsilon2
    
    # Here we append the output of the echo command to the logfile using "tee".
    # The -a flag means "append".
    echo "Submitted job batch #$i with Job ID $jobid, EPSILON1=$EPSILON1, EPSILON2=$EPSILON2" | tee -a $logfile
	
	./slip_sim_cluster -t 100000 -e $epsilon1 -r $epsilon2 --jobid $SLURM_JOB_ID #-o 1

	echo "Job batch #$i (Job ID: $SLURM_JOB_ID) with EPSILON1=$EPSILON1, EPSILON2=$EPSILON2 completed at $(date)" | tee -a $logfile
	
	i=$(($i + 1))

done < "$input_file"  # Read the input file line by line

echo "All jobs submitted. Monitoring in the background. Check $logfile for updates."

# Summary of batch job info saved to a file (Optional).
qstat -f $SLURM_JOB_ID > $SLURM_JOB_NAME-$SLURM_JOB_ID.info-summary
