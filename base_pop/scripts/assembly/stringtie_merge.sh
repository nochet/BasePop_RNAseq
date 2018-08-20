#!/bin/bash
#-------------------------------------------------------------------------------
#  SBATCH CONFIG
#-------------------------------------------------------------------------------
## resources
#SBATCH --partition Lewis
#SBATCH --nodes=1
#SBATCH --ntasks=1  # used for MP#SBATCH -e error_%A_%a.err # Standard errorI codes, otherwise leave at '1'
#SBATCH --cpus-per-task=12  # cores per task
#SBATCH --mem-per-cpu=8G  # memory per core (default is 1GB/core)
##SBATCH --array=1-162
#SBATCH --time 0-02:00  # days-hours:minutes
#SBATCH --qos=normal
#SBATCH --account=kinglab  # investors will replace this with their account name
#
## labels and outputs
# labels and outputs
#
#SBATCH --job-name=StringTie_merge_step4
#
#SBATCH -o test_%A_%a.out # Standard output
#SBATCH -e error_%A_%a.err # Standard error


#
## notifications
#SBATCH --mail-user=ngomae@missouri.edu  # email address for notifications
#SBATCH --mail-type=END,FAIL  # which type of notifications to send
#-------------------------------------------------------------------------------

echo "### Starting at: $(date) ###"

module load stringtie/stringtie-1.3.3b

COMMANDA=`head -n ${SLURM_ARRAY_TASK_ID} mergelist.txt | tail -n 1`
eval $COMMANDA

echo "### Ending at: $(date) ###"

