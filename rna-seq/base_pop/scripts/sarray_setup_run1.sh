#!/bin/bash
#-------------------------------------------------------------------------------
#  SBATCH CONFIG
#-------------------------------------------------------------------------------
## resources
#SBATCH --partition BioCompute
#SBATCH --nodes=1
#SBATCH --ntasks=1  # used for MP#SBATCH -e error_%A_%a.err # Standard errorI codes, otherwise leave at '1'
##SBATCH --cpus-per-task=12  # cores per task
#SBATCH --mem-per-cpu=8G  # memory per core (default is 1GB/core)
#SBATCH --time 0-02:00  # days-hours:minutes
#SBATCH --qos=normal
#SBATCH --account=kinglab  # investors will replace this with their account name
#
## labels and outputs
#SBATCH --job-name=rna1_enoch
#SBATCH --output=results-%j.out  # %j is the unique jobID
#
#SBATCH -o test_%A_%a.out # Standard output
#SBATCH -e error_%A_%a.err # Standard error

## notifications
#SBATCH --mail-user=ngomae@missouri.edu  # email address for notifications
#SBATCH --mail-type=END,FAIL  # which type of notifications to send
#-------------------------------------------------------------------------------

echo "### Starting at: $(date) ###"


# load packages
module load hisat2/hisat2-2.0.5

COMMANDA=`head -n ${SLURM_ARRAY_TASK_ID} scripts/Align_cmd_run1.txt | tail -n 1`
$COMMANDA


echo "### Ending at: $(date) ###"
 
