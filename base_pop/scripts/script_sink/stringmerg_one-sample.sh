#!/bin/bash
#-------------------------------------------------------------------------------
#  SBATCH CONFIG
#-------------------------------------------------------------------------------
## resources
#SBATCH --partition BioCompute
#SBATCH --nodes=1
#SBATCH --ntasks=1  # used for MP#SBATCH -e error_%A_%a.err # Standard errorI codes, otherwise leave at '1'
#SBATCH --cpus-per-task=12  # cores per task
#SBATCH --mem-per-cpu=8G  # memory per core (default is 1GB/core)
#SBATCH --time 0-02:00  # days-hours:minutes
#SBATCH --qos=normal
#SBATCH --account=kinglab  # investors will replace this with their account name
#
## labels and outputs
#SBATCH --job-name=stringmerg_one-sample
#SBATCH --output=results-%j.out  # %j is the unique jobID
#
#SBATCH -o stringmerg_%a.out # Standard output
#SBATCH -e error_%a.err # Standard error

## notifications
#SBATCH --mail-user=ngomae@missouri.edu  # email address for notifications
#SBATCH --mail-type=END,FAIL  # which type of notifications to send
#-------------------------------------------------------------------------------

echo "### Starting at: $(date) ###"

# load modules then display what we have

## module avail -t stringtie shows several including stringtie/stringtie-1.3.3b

module load stringtie/stringtie-1.3.3b

stringtie --merge -p 8 -G genes/dmel-all-r6.18.gtf -o stringtie_merged.gtf processed/merglist.txt

echo "### Ending at: $(date) ###"

