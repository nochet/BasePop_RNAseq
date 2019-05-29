
---
title: "RNAseq Base Population Project Record & Notes"
author: "Enoch Ng'oma (EN), and Elizabeth King (EK)"
date: "7/9/2018"
---


## Project Record

- Start at the top of the Methods section.
- When you get to a data analysis step, identify the script that produced it and put an entry in to the project record. You can use subsections- try to use same as in methods or results.
- Then go through results. Make sure each result has a source entry for the script that made it, including visualizations. 
- Also, do all of these run? Are the paths correct?

e.g. 

## RNA-seq processing
Step 1:....

Step 2:....




Genome-wide expression (RNA-seq) of base population

A. download reference files (EN)
gtf files downloaded from Index of ftp://ftp.flybase.org/genomes/dmel/current/

B. Sequence assembly (ref. Pertea et al 2016)
`Set_up_arrays.Rmd` - step-by-step protocol as follows:



STEP 1: Hisat2 (Align reads to reference)

Test: align one sample to genome:
hisat2 --dta -q -x indexes/bdgp6_tran/genome_tran -U samples/run1/HS-6_O_S53_R1_001.fastq.gzHS-6_O_S53_R1_001.fastq.gz -S HS6O_1.sam >temp.txt 2>error.txt &

Do Step 1 Align in `/base_pop/scripts/assembly_longProtocol/Set_up_arrays.Rmd`
Output files: Align_cmd_run1.txt, Align_cmd_run2.txt and Align_cmd_run3.txt 
Copy output .txt files to Lewis Cluster

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


STEP 2b: Samtools merge (combine transcripts from all runs)

Do Samtools merge step in `Set_up_arrays.Rmd` by running `???_merge.sh` check lewis
Output: `S02b_samtools_merge.txt`
Next, run `samtools_merge.sh` to merge


STEP 4: StringTie (Assemble and quantify expressed genes and transcripts)

Do StringTie (Assemble) in `Set_up_arrays.Rmd` by running `???.sh` see lewis
Output: `S04_stringtie_assemble.txt`
Next, run `stringtie_assemble.sh`
Output: ?????

Run perl script to extract more gene ids from StringTie output

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


STEP8A: BallGown (transcript level differential expression based on transformed transcript or gene counts)

STEP8B: DESeq (gene-level differential expression based on read count)
Run `prepDESeq.Rmd` by following instructions in `/scripts/DESeq_scripts/prepDEpy_instructions.txt`

STEP9A: Control for batch effects using SVAseq package
- `dataPrep_batch.Rmd` - data prep for sva
	- input1: `/processed/describe_samples_batch.csv` 
	- input2: `/processed/results/ballG_all_results/bg_ballG_all_results.Rda`
	- output1: `/processed/results/ballG_all_results/log_gfpkm_all.Rda` log-transformed data for sva
	- output2: `/processed/results/ballG_all_results/nolog_gfpkm_all.Rda` untransformed data for svaseq
- `DiffExpr_batch.Rmd` - identify and remove batch effects
	- input: `/processed/results/ballG_all_results/nolog_gfpkm_all.Rda`
	- output: `svaseq.dat` - an object containing batch-corrected expression data from the sva package
	
STEP9B: 

## Project Notes

### General Notes & links

 Get more gene_ids from StringTie object with python: https://www.biostars.org/p/261128/

 Ref: http://ccb.jhu.edu/software/stringtie/index.shtml?t=manual#de
 Ref p.92-93 Haddok & Dunn book
 Detailed instructions are in /scripts/DESeq_scripts/prepDEpy_instructions.txt

 To split a .csv column with multiple id entries, follow: http://help.snapizzi.com/csv-files/split-csv-data-into-different-columns

DESeq vs Ballgown, see https://support.bioconductor.org/p/107011/

Precomputed FlyBase files:
http://flybase.org/cgi-bin/get_static_page.pl?file=bulkdata7.html&title=Current%20Release




2018-12-27 (EN)
Wanted to compare results with a different approach: ENRICHMENTBROWSER
`/scripts/GOscripts/enrichBrowse_GO.Rmd`
Stuck on achieving a `SummarizedExperiment` class matrix of results
I have a `DESeqResults` class object

2018-12-22 (EN)
Gene set enrichment analysis using GAGE
`/scripts/GOscripts/gsea_GO.Rmd`
	Load `load("../../processed/DESEQ/GO/lrt.treatment.Rda")`
	Need libraries: gage, pathview and org.Dm.eg.db
	Followed this tutorial alongside all 3 vignettes: http://www.gettinggeneticsdone.com/2015/12/tutorial-rna-seq-differential.html
		Mapped FlyBase CG SYMBOLs in our results to fly annotation database `org.Dm.eg.db`
		Converted CG SYMBOLs to ENTREZ
		Performed pathway and GO analyses

2018-11-30 (EGK)
Computed fold changes with shrinkage (method= normal). Made plots of fold changes in DR and HS relative to control for body, ovary & head. P-values are the adjusted p-values for the OVERALL effect of treatment, so could be driven by one diet or tissue. Script is foldchange_DEseq.Rmd, fold change data is in all_fc_dat.rda and plot is FC_all.pdf


2018-11-07 (EN)
Reduced number of scripts to 3, renamed as follows:
	`prep_DESeq.Rmd` renamed `batch_DESeq.Rmd` - runs SVA
	`bg_shortP_viz.Rmd` renamed `short_protc_viz.Rmd` - makes visualizations in ballgown and DESeq2
	`DESeq_DE.Rmd` renamed `DESeq_DExpr.Rmd` - fits models for differential expression in DESeq2 based expression values from the short protocol of stringtie
	

2018-11-02 (EN)
PCA plots in DESeq2: `DESeq_viz.Rmd`
Input1 `/processed/DESEQ/sampleDat_with_SV.csv` (`made in prep_DESe2.Rmd`)
Input2: `/processed/DESEQ/normExpr_countData.csv` (sva batch-controlled short protocol expression `made in prep_DESe2.Rmd`)
Run `deseq()` on model without- and model with interaction
Transform DESeq results using 	1) Variance Stabilization Transformation (vst),
								2) Regularized Log Transformation (rlog)
Plot 2 PCA graphs for PC1/PC2, one for for each transformation; both for the interaction because PCA plots for both DESeq models not different 

2018-10-26 (EN)
Making PCA for DESeq data
	Get fpkm values with fpkm(). An error occurs. Solution is to use tximport package to normalize count values for gene size https://support.bioconductor.org/p/83607/
	Starts from "stringtie -eB -G transcripts.gff <source_file.bam>"
	cd to /processed/S03_short_ballG/ballgown: 
		find * -type d > ../sampNames.txt
		write $ in the search window i.e. matches end of line
		write _ballgown.gtf in the replace window
		write ^ in search window i.e. matches start of line
		write ../../processed/S03_short_ballG/ballgown

2018-10-25 (EN)
DEGs for treatment effect - 2475 genes
DEGs for interaction effect - 977 genes
Identified 3 DEGs under trans-eQTL (step and Ilp5) peaks. In script `genes_in_QTL.Rmd`, see object `distt`
	- Ahcy (step QTL) is found in both treatment and interaction lists. Gene summary file associates it with: viable, aging, fecundity. Involved in NAD processes
	- Cda4 (Ilp5 QTL) occurs only in treatment list. Involved in chitin and carbohydrate metabolism via NADH activity
	- CG40486 in treatment list
slif was not DE in both treatment and interaction effects.

### 2018-10-21 (EN)
Scanned Stanley et al 2017 trans-eQTL and lifespan QTL for DEGs - see R object `qtl_degs`
Split `prep_DESeq.Rmd` script into:
	`prep_DESeq.Rmd`: runs batch control of Stringtie output using SVA
	`BaseP_DESeq.Rmd`: Differential expression models with DESeq2
	`BaseP_DESeq_Visuals.Rmd`: plots to visualize expression data

### 2018-10-20
`genes_in_QTL.Rmd` 
Code to scan past QTL peaks (Stanley et al 2017) for DEGs (this study)
Reads and pulses `"../../processed/gene_map_table_fb_2018_05.tsv"` downloaded from ftp://ftp.flybase.net/releases/FB2018_05/precomputed_files/genes/. 
A function, `DEunderPeak()` scans under QTL peaks.
A prepared data frame is saved to `"../../processed/DESEQ/DEG_QTL/gene_map_table.csv")`.
Steps:
	- merge a dataframe of DEGs (`DEG_treat`) and the gene map table (`gmtable`) on FBgn
	- read in data with QTL peaks `../../processed/DESEQ/DEG_QTL/iis_table1.txt`. This excludes peak that span the centromere
	- scan through trans eQTL and lifespan peaks for DEGs

### 2018-10-10 (EN)
Download a precomputed gene annotation file for use in GO analysis
Fbgn annotation file "Genes:fbgn_annotation_ID_fb_2018_05.tsv"
Index of ftp://ftp.flybase.net/releases/FB2018_05/precomputed_files/genes/

### 2018-10-09
Convert FBgn symbols to name and annotation symbols for use with GO analysis packages
Go to FlyBase ID Converter tool: http://flybase.org/convert/id
Paste list of FBgn.. and choose 'Validate and Convert' to 'Genes'
On result page choose 'BatchDownload'
Select Symbol, Name, Annotation Symbol, FlyBase ID. Push 'Get Field Data'
Press 'Download as a file'

### 2018-10-07 (EN)
DAVID web tool for Gene Ontology Analysis. I used the Functional Annotation Tool to get GO terms in all 3 categories (BP, MF and CC). Also obtained annotation chart and table
969 genes sig differentiated (padj<0.05, LRT) from the interaction model, `dds.02`

`GO_deseq.Rmd`
GOplot to visualize: https://cran.r-project.org/web/packages/GOplot/vignettes/GOplot_vignette.html 


### 2018-10-02
2) Use GOrilla for functional analysis of each list of DEGs a) as a whole and b) up- down- separately
Discontinued.

### 2018-09-04 (EN)
Run a model in DESeq2 to test effect of tissue x treatment interaction on normalized batch-controlled gene-wise data:
dds.02 <- DESeqDataSetFromMatrix(countData = countdata, 
                                   colData = phenDat.sva, 
                                   design = ~ SV1 + batch + tissue*treatment)

The following pairwise comparisions for treatment effect:
rr1 <- results(dds_deseq.02, contrast=c('tissue', "H","B"), alpha = 0.05)
rr2 <- results(dds_deseq.02, contrast=c('tissue', "O","B"), alpha = 0.05)
rr3 <- results(dds_deseq.02, contrast=c('tissue', "O","H"), alpha = 0.05)

and,

tissue effect:
ff1 <- results(dds_deseq.02, contrast=c('treatment', "C","DR"), alpha = 0.05)
ff2 <- results(dds_deseq.02, contrast=c('treatment', "C","HS"), alpha = 0.05)
ff3 <- results(dds_deseq.02, contrast=c('treatment', "DR","HS"), alpha = 0.05)

resulted in tables of DEGs.

1) I compiled a list of 26 genes previously identified in QTL studies of diet treatments (Stanley et al, 2017): `/processed/DESEQ/DEG_QTL/iis_genelist.csv"`. I identify any of these which are also significantly expressed in this study: files written to `/processed/DESEQ/DEG_QTL/`

Genes identified are collected in one file: `/processed/DESEQ/DEG_QTL/qtl_genes_all_FlyBase.csv`. Also present in this file is FlyBase genome and functional data for each gene from a batch search. 

### 2018-08-30 (EGK)

Developed compare_StringTie_method script. Compares: 1) gene count table generated from stringtie run that includes estimating novel transcripts, 2) gene count table generated from stringtie run that includes estimating novel transcripts followed by perl script to get back FBgn names, 3) gene count table ignoring novel transcripts

General comparison:

Type 1 has very few FBgn names. When comparing type 2 and 3 (merging on gene names), those with differences are almost all a MSTRG that is split up in the perl script to assign gene names. So some reads are assigned to a gene and some just get an MSTRG number (both share the same string tie MSTRG number). Another set are genes with more than one assigned with the perl script. It makes sense these do not agree. Those that are left ~800 I am not sure why they don't agree. For one I looked at, the pattern matches what modEncode reports on flybase (low ovary expression) for type 3, while type 2 and 1 are both way higher for all tissues including ovaries. In conclusion, I feel good about using type 3 for GENE level analysis. 

Only 200 gene names in type 3 are not represented in type 2. Not sure why for these.

### 2018-08-27 (EN)

Implementing an alternative workflow. A visual of the short protocol is here: http://ccb.jhu.edu/software/stringtie/index.shtml?t=manual#de

Run stringtie -eB on the output of samtools -sort ie. (step 2 of Pertea et al, 2016)

Gene count tables are written to /processed/S03_short_ballG/ `gene_count_matrix` and `transcript_count_matrix`

### 2018-08-27 (EGK)

Attempting to generate novel transcripts and get an accurate gene level analysis is complex. We plan to perform separate pipelines for a transcript level and gene level analysis. 

Gene level analysis: we will use the reference annotation rather than performing the merge step in StringTie followed by DEseq

Transcript level analysis: we will use the StringTie-Ballgown pipeline

Exploration of the gene count table generated from the data using the merge file before and after assigning a reference gene gave different results. Some MSTRGs are split up into those that are assigned to a gene and those that are not. 

### 2018-08-20 (EN)

- Get more gene ids tranferred from StringTie by running a perl script `mstrg_prep.pl` explained here:
https://github.com/gpertea/stringtie/issues/179
and documented here:
https://gist.github.com/gpertea/b83f1b32435e166afa92a2d388527f4b
- Command: `perl mstrg_prep.pl stringtie_merged.gtf > stringtie_merged_pl.gtf`
- Then, re-run step 6 (Pertea et al, 2016) to compute transcript abundances:  sbatch `stringtie_abundances.sh` 
- Results stored in `ballG_pl`
- Problem: duplicated MSTRG and Fbgn ids. Perl script appears to assign transcript if not sure, to the same MSTRG id creating duplicates.

### 2018-07-30 (EN)

- Implemented DESeq2 on svaseq-prepped read count data
- 22 rows did not converged when running `ddsDE.ttSS <- DESeq(dds.ttSS)`
- Dealing with autocorrelation: see https://support.bioconductor.org/p/65091/ and tried:
- 1) Filter out normalized count of at least 10 reads in two or more samples. 
	- `filter <- rowSums(nc >= 10) >= 2`
	- 8 rows did not converge
- 2) Try increasing the 2 sample requirement to 3 or 4
	- `filter <- rowSums(nc >= 10) >= 3`
	- 6 rows did not converge
- 3) omit non-converging rows from the results step
	- `ddsDE_Clean <- ddsDE_filt10[which(mcols(ddsDE_filt10)$betaConv),]`
	- didn't run, but also doesn't feel a good approach!
- 4) increase the maximum iterations (maxit) instead of running DESeq()
	- `filter <- rowSums(nc >= 10) >= 3` plus`nbinomWaldTest(dds, maxit=500)`
- 5) `filter <- rowSums(nc >= 5) >= 2` plus `nbinomWaldTest(maxit=500)`
	- 9 rows don't converge

### 2018-07-24 (EGK)

- Looked through DEseq code to check
- Implemented svaseq to check for batch effects
    - Will follow this: https://support.bioconductor.org/p/80755/
    - Needed to remove all genes with sample counts that were all zeros. (gave this error: Error in density.default(x, adjust = adj) : 'x' contains missing values)
    - sva identified 2 svas
    - clearly 4 groupings when plotting svas- unknown sources

### 2018-07-20 (EN)

- Performed first DESeq model: `design = ~ treatment + tissue`. Need to 1) do batch analysis and 2) include an interaction term `treatment*tissue`

### 2018-07-19 (EN)

- Sub-structured the `scripts` directory into `Ballgown_scripts` and `DESeq_scripts`
- Renamed `/processed/ballG_all/` as ` /processed/ballgown`. The directory structure created bt StringTie is required `prepDE.py` to compile counts (see http://ccb.jhu.edu/software/stringtie/index.shtml?t=manual#de)

### 2018-07-17 (EN)

- Started a new, hopefully streamlined `DiffExpr_batch.Rmd` script. Script with all models tried up to now is pushed to `script_sink/` an `DiffExpr_batch_old.Rmd`

### 2018-07-11 (EN)

- Move all ballgown stuff from `DiffExpr_batch.Rmd` to `DiffExpr.Rmd`
- clean redundant code (i.e. same code in multiple scripts)
- change column `group` in `phenDat` to `batch`

### 2018-07-09 

Plan for analysis: (EN)

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


















