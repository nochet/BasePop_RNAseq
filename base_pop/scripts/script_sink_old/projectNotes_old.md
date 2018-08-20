Genome-wide expression (RNA-seq) of base population

1. download reference files (EN)
gtf files downloaded from Index of ftp://ftp.flybase.org/genomes/dmel/current/

Setup alignment on nivalis server
test: align one sample to genome (ref. Pertea et al 2016):
hisat2 --dta -q -x indexes/bdgp6_tran/genome_tran -U samples/run1/HS-6_O_S53_R1_001.fastq.gzHS-6_O_S53_R1_001.fastq.gz -S HS6O_1.sam >temp.txt 2>error.txt &

Interactively: Create command list for each of 54 alignments in each run:
cd to scripts/ folder when running code
```{r}
# code reads file names from Libby's directory on nivalis and creates all 162 alignment commands, saving them to three files named by run (1-3)
# run code in R studio in a web browser - not commandline R
runs<-c(1,2,3)
for(run in runs)
{
ffs <- list.files(paste("/home/kingeg/Projects/RNA_SEQ_BasePop/Data/run",run,"/",sep=""), pattern="fastq.gz$")

for(ii in 1:length(ffs))
{
st.s <- "hisat2 --dta -q -x indexes/bdgp6_tran/genome_tran -U " 

st.p <- paste("samples/run",run,"/",sep="")

outf <- paste("processed/",strsplit(ffs[ii], ".", fixed=TRUE)[[1]][1],"_run",run,sep="")

cat(paste(st.s, st.p, ffs[ii], " -S ", outf,".sam",sep=""),"\n", file=paste("scripts/Align_cmd_run",run,".txt",sep=""),append=TRUE)

}

}

# Note: code works if I cd to scripts directory; if not there, change paths accordingly

# copy the output files: Align_cmd_run1.txt, Align_cmd_run2.txt and Align_cmd_run3.txt to Lewis Cluster
```


### Prepare one sbatch file (.sh) for each run - script: 'alig_setup.sh'
# example here is for run1; change names and parameters accordingly

!/bin/bash
#-------------------------------------------------------------------------------
#  SBATCH CONFIG
#-------------------------------------------------------------------------------
## resources
#SBATCH --partition BioCompute
#SBATCH --nodes=1
#SBATCH --ntasks=1  # used for MP#SBATCH -e error_%A_%a.err # Standard errorI codes, otherwise leav$
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
## this command reads commands from 'Align_cmd_run1.txt' and runs alignment

echo "### Ending at: $(date) ###"


# run test align of two samples to see if it works
sbatch --array=1-2 scripts/sarray_setup_run1.sh

# Now, run for all samples in a run1
sbatch --array=1-54 scripts/sarray_setup_run1.sh

# can run the three sbatch files at once - the scheduler takes care of resources



## This sbatch file 'alig-sample.sh' aligns one sample interactively on Lewis

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
#SBATCH --job-name=alig_test
#
#SBATCH -o test_%a.out # Standard output
#SBATCH -e error_%a.err # Standard error

## notifications
#SBATCH --mail-user=ngomae@missouri.edu  # email address for notifications
#SBATCH --mail-type=END,FAIL  # which type of notifications to send
#-------------------------------------------------------------------------------

echo "### Starting at: $(date) ###"

# load modules then display what we have
module load hisat2/hisat2-2.0.5
# module avail hisat2 shows hisat2/hisat2-2.0.5

hisat2 --dta -q -x indexes/bdgp6_tran/genome_tran -U samples/run3/HS-3_O_S45_R1_001.fastq.gzHS-3_O_S45_R1_001.fastq.gz -S processed/HS-3_O_S45_R1_001_run3.sam

echo "### Ending at: $(date) ###"
###################


sarray_setup_run1.sh  # for 54 samples in run1
sarray_setup_run2.sh  # 54 samples in run2
sarray_setup_run2.sh  # 54 samples in run3



### Sort and convert the SAM files to BAM
# example from hisat2 paper
$ samtools sort -@ 8 -o ERR188044_chrX.bam ERR188044_chrX.sam

```{r}

#samtools sort -@ 8 -o C-3_O_S9_R1_001_run2.bam C-3_O_S9_R1_001_run2.sam 
# samtools sort -o processed/HS-3_O_S45_R1_001_run3.bam processed/HS-3_O_S45_R1_001_run3.sam

ll<-list.files("processed/", pattern="sam$")

for(ii in ll)
{
   cc<-strsplit(ii,".",fixed=TRUE)[[1]][1]
   cat(paste("samtools merge -r/",cc,".bam ", "processed/", ii,sep=""), "\n", file=paste("scripts/Align_cmd_run",run,".txt",sep=""),append=TRUE)
   
  
}

```


# One sample: 'samtools_one-sample.sh'

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

# All 162 .sam to .bam
# samtools_setup.sh runs R-script samtools_all.R; output is in /scripts/samtools_all.txt
# sacct -j 4212835 --format=JobID,JobName,MaxRSS,Elapsed # check time for 1 sample (i.e. 57s)

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
#SBATCH --job-name=rna_enoch
#
#SBATCH -o test_%A_%a.out # Standard output
#SBATCH -e error_%A_%a.err # Standard error

## notifications
#SBATCH --mail-user=ngomae@missouri.edu  # email address for notifications
#SBATCH --mail-type=END,FAIL  # which type of notifications to send
#-------------------------------------------------------------------------------

echo "### Starting at: $(date) ###"


# load packages
module load samtools/samtools-1.3.1

COMMANDA=`head -n ${SLURM_ARRAY_TASK_ID} samtools_all.txt | tail -n 1`
$COMMANDA

echo "### Ending at: $(date) ###"

#### sbatch --array=1-162 samtools_all.sh ####



#### Assemble transcripts for one sample: sbatch stringt_one-sample.sh
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
#SBATCH --job-name=stringt_one-sample
#SBATCH --output=results-%j.out  # %j is the unique jobID
#
#SBATCH -o stringt_%a.out # Standard output
#SBATCH -e error_%a.err # Standard error

## notifications
#SBATCH --mail-user=ngomae@missouri.edu  # email address for notifications
#SBATCH --mail-type=END,FAIL  # which type of notifications to send
#-------------------------------------------------------------------------------

echo "### Starting at: $(date) ###"

# load modules then display what we have

## module avail -t stringtie shows several including stringtie/stringtie-1.3.3b

module load stringtie/stringtie-1.3.3b

stringtie -p 8 -G genes/dmel-all-r6.18.gtf -o processed/HS-3_O_S45_R1_001_run3.gtf -l HS-3_O_S45_R1_001_run3 processed/HS-3_O_S45_R1_001_$

echo "### Ending at: $(date) ###"

### For all samples

2018-06-19 (EN)
# Tested for batch effects using SVA package. The "be" and "leek" methods applied ns.v() produced 11 and 0 SVs, respectively.Wondering why the difference we read:
  # Leek_2014_ svaseq: removing batch effects and other unwanted noise from sequencing data, 
  # Leek_2011_Asymptotic Conditional Singular Value Decomposition for High-Dimensional Genomic Data,
  # https://support.bioconductor.org/p/97469/,
  # Jaffe et al_2015_ Practical impacts of genomic data “cleaning” on biological discovery using surrogate variable analysis
  # Found that methods can produce diferent results

2018-06-20 (EN & EGK)
# Decided: 1) check that known batches in our case have effect on expression (a. include batch in PCA, b. test with sva or linear model) 2) run svaseq without known bathes

















