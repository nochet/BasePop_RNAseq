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
#SBATCH --job-name=alig_setup
#
#SBATCH -o setup_%A_%a.out # Standard output
#SBATCH -e error_%A_%a.err # Standard error

## notifications
#SBATCH --mail-user=ngomae@missouri.edu  # email address for notifications
#SBATCH --mail-type=END,FAIL  # which type of notifications to send
#-------------------------------------------------------------------------------

echo "### Starting at: $(date) ###"

module load R/R-3.3.3
R
# code adapted from egk
runs<-c(1,2,3)
for(run in runs)
{
ffs <- list.files(paste("../samples/run",run,"/",sep=""), pattern="fastq.gz$")

for(ii in 1:length(ffs))
{
st.s <- "hisat2 --dta -q -x indexes/bdgp6_tran/genome_tran -U " 

st.p <- paste("samples/run",run,"/",sep="")

outf <- paste("processed/",strsplit(ffs[ii], ".", fixed=TRUE)[[1]][1],"_run",run,sep="")

cat(paste(st.s, st.p, ffs[ii], " -S ", outf,".sam >temp",ii,".txt 2>error",ii,".txt &",sep=""),"\n", file=paste("Align_cmd_run",run,".txt",sep=""),append=TRUE)

}

}
# Note: code works if I cd to scripts directory; if not there, change paths accordingly

echo "### Ending at: $(date) ###"
