#!/bin/bash

logfile="flare_batch_submission.log"
input_file="input.txt"  # Change this to your input file with numerical values
i=1
# Note that $(something) calls the 'something' command and replaces $(something) with the output.
# Note that $something replaces $something with the value of the 'something' variable
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


     # Submit the job array and capture its job ID
    jobid=$(sbatch slip_sim_bootstraps.sbatch | awk '{print $4}')
    # Here we append the output of the echo command to the logfile using "tee".
    # The -a flag means "append".
    echo "Submitted job batch #$i with Job ID $jobid, EPSILON1=$EPSILON1, EPSILON2=$EPSILON2" | tee -a $logfile

    # Start a background process to monitor job completion
    (
        while squeue -t pending,running -j $jobid 2>/dev/null | grep -q $jobid; do
            sleep 30  # Check every 30 seconds
        done
        echo "Job batch #$i (Job ID: $jobid) with EPSILON1=$EPSILON1, EPSILON2=$EPSILON2 completed at $(date)" | tee -a $logfile
    ) &  # Run in the background
	
	i=$(($i + 1))

done < "$input_file"  # Read the input file line by line

echo "All jobs submitted. Monitoring in the background. Check $logfile for updates."