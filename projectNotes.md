
---
title: "RNAseq Base Population Project Record & Notes"
author: "Enoch Ng'oma (EN), and Elizabeth King (EK)"
date: "7/9/2018"
---


## Project Record



Genome-wide expression (RNA-seq) of base population

A. download reference files (EN)
gtf files downloaded from Index of ftp://ftp.flybase.org/genomes/dmel/current/

B. Sequence assembly (ref. Pertea et al 2016)
`Set_up_arrays.Rmd` - step-by-step protocol as follows:


STEP 1: Hisat2 (Align reads to reference)

test: align one sample to genome:
hisat2 --dta -q -x indexes/bdgp6_tran/genome_tran -U samples/run1/HS-6_O_S53_R1_001.fastq.gzHS-6_O_S53_R1_001.fastq.gz -S HS6O_1.sam >temp.txt 2>error.txt &

Do Step 1 Align in `Set_up_arrays.Rmd`
Output files: Align_cmd_run1.txt, Align_cmd_run2.txt and Align_cmd_run3.txt 
Copy output files to Lewis Cluster

Commands at the prompt:
first run test align of two samples to see if it works:
sbatch --array=1-2 scripts/sarray_setup_run1.sh

second, run for all samples in a run1:
sbatch --array=1-54 scripts/sarray_setup_run1.sh

sbatch scripts:
`sarray_setup_run1.sh`  # 54 samples in run1
`sarray_setup_run2.sh`  # 54 samples in run2
`sarray_setup_run2.sh`  # 54 samples in run3

Output files: ???.sam (SEE Lewis directory), path???

Note: can run the three sbatch files at once - the scheduler takes care of resources


STEP 2: Samtools (Sort and convert the SAM files to BAM)

example from hisat2 paper:
$ samtools sort -@ 8 -o ERR188044_chrX.bam ERR188044_chrX.sam

Do Samtools sort in `Set_up_arrays.Rmd` by running `samtools_setup.sh`  
Output: `S02_samtools_all.txt`
Next, run `samtools_all.sh`
Output: ???dot bams

Command at the prompt:
sbatch --array=1-162 samtools_all.sh


STEP 3: Samtools merge (combine transcripts from all runs)

Do Samtools merge step in `Set_up_arrays.Rmd` by running `???_merge.sh` check lewis
Output: `S03_samtools_merge.txt`
Next, run `samtools_merge.sh` to merge


STEP 4: StringTie (Assemble and quantify expressed genes and transcripts)

Do StringTie (Assemble) in `Set_up_arrays.Rmd` by running `???.sh` see lewis
Output: `S04_stringtie_assemble.txt`
Next, run `stringtie_assemble.sh`
Output: ?????

STEP 5: gffcompare (Compare transcripts to reference genome - optional)
step skipped

STEP 6: StringTie (Transcript abundances and table of counts)

Do StringTie transcript abundances in `Set_up_arrays.Rmd` by running `???.sh` see lewis
Output: `S06_abundances.txt`
Next, run `stringtie_abundances.sh`
Output: ????

STEP 7: Create a csv containing sample ids

Do "Create a csv containing sample ids" in `Set_up_arrays.Rmd`
Output: `describe_samples.csv`


STEP8: BallGown (Differential expression)

STEP9: Control for batch effects using SVAseq package
- `dataPrep_batch.Rmd` - data prep for sva
	- input1: `/processed/describe_samples_batch.csv` 
	- input2: `/processed/results/ballG_all_results/bg_ballG_all_results.Rda`
	- output1: `/processed/results/ballG_all_results/log_gfpkm_all.Rda` log-transformed data for sva
	- output2: `/processed/results/ballG_all_results/nolog_gfpkm_all.Rda` untransformed data for svaseq
- `DiffExpr_batch.Rmd` - identify and remove batch effects
	- input: `/processed/results/ballG_all_results/nolog_gfpkm_all.Rda`
	- output: `svaseq.dat` - an object containing batch-corrected expression data from the sva package
	


## Project Notes 

### 2018-07-11 (EN)
- Move all ballgown stuff from `DiffExpr_batch.Rmd` to `DiffExpr.Rmd`
- clean redundant code (i.e. same code in multiple scripts)
- change column `group` in `phenDat` to `batch`

### 2018-07-09 

Plan for analysis:

Compare ballgown and sva differential expression analyses
 
- Ballgown done in `DiffExpr.Rmd`
- Prep data for sva - `dataPrep_batch.Rmd`
- Do sva - `DiffExpr_batch.Rmd`
- Compare sva & ballgown

### 2018-06-20 (EN & EGK)
- Decided: 1) check that known batches in our case have effect on expression (a. include batch in PCA, b. test with sva or linear model) 2) run svaseq without known bathes


### 2018-06-19 (EN)
- Tested for batch effects using SVA package. The "be" and "leek" methods applied ns.v() produced 11 and 0 SVs, respectively.Wondering why the difference we read:
  - Leek_2014_ svaseq: removing batch effects and other unwanted noise from sequencing data, 
  - Leek_2011_Asymptotic Conditional Singular Value Decomposition for High-Dimensional Genomic Data,
  - https://support.bioconductor.org/p/97469/,
  - Jaffe et al_2015_ Practical impacts of genomic data “cleaning” on biological discovery using surrogate variable analysis
  - Found that methods can produce diferent results


















