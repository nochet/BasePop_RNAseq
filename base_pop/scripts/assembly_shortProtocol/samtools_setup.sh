#!/bin/bash
#-------------------------------------------------------------------------------
#  SBATCH CONFIG
#-------------------------------------------------------------------------------
## resources
#SBATCH --partition BioCompute
#SBATCH --nodes=1
#SBATCH --ntasks=1  # used for MP#SBATCH -e error_%A_%a.err # Standard errorI codes, otherwise leave at '1'
#SBATCH --cpus-per-task=1  # cores per task
#SBATCH --mem=100  # memory per core (default is 1GB/core)
#SBATCH --time 0-00:10  # days-hours:minutes
#SBATCH --qos=normal
#SBATCH --account=kinglab  # investors will replace this with their account name
#
## labels and outputs
#SBATCH --job-name=sam_setup
#
#SBATCH -o setup_%A_%a.out # Standard output
#SBATCH -e error_%A_%a.err # Standard error

## notifications
#SBATCH --mail-user=ngomae@missouri.edu  # email address for notifications
#SBATCH --mail-type=END,FAIL  # which type of notifications to send
#-------------------------------------------------------------------------------

echo "### Starting at: $(date) ###"
module load R/R-3.3.3

R  --no-save -f samtools_all.R


echo "### Ending at: $(date) ###"

