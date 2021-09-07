#!/bin/csh
#SBATCH -J CapSol
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --account=PSS0003

# A Basic Gaussian Job for the OSC Owens Cluster
#   https://www.osc.edu/resources/available_software/software_list/gaussian

# SLURM_SUBMIT_DIR refers to the directory from which the job was submitted.

# Emit verbose details as the job script executes.
set echo
set verbose
# Emit verbose details on the job's queuing.
squeue -j $SLURM_JOB_ID


# The following lines set up the Capsol environment



#
# Move to the directory where the job was submitted
#
cd $TMPDIR
cp -p "$SLURM_SUBMIT_DIR/RR_originalTirmzi-0003.in" capsol.in
#
# Run Orca
#
/fs/ess/PSS0003/bin/capsol > "$SLURM_SUBMIT_DIR/RR_originalTirmzi-0003.in.out"
ls -al
#
# Now, copy data (or move) back once the simulation has completed
cp * $SLURM_SUBMIT_DIR