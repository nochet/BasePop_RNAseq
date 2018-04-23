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
#SBATCH --time 0-00:40  # days-hours:minutes
#SBATCH --qos=normal
#SBATCH --account=kinglab  # investors will replace this with their account name
#
## labels and outputs
#SBATCH --job-name=samtools_one-sample
#SBATCH --output=results-%j.out  # %j is the unique jobID
#
#SBATCH -o test_%a.out # Standard output
#SBATCH -e error_%a.err # Standard error

## notifications
#SBATCH --mail-user=ngomae@missouri.edu  # email address for notifications
#SBATCH --mail-type=END,FAIL  # which type of notifications to send
#-------------------------------------------------------------------------------

echo "### Starting at: $(date) ###"

# load modules then display what we have

## module avail -t samtools shows several including samtools/samtools-1.3.1 

module load samtools/samtools-1.3.1

samtools sort -@ 8 -o ../processed/HS-3_O_S45_R1_001_run3.bam ../processed/HS-3_O_S45_R1_001_run3.sam

echo "### Ending at: $(date) ###"
