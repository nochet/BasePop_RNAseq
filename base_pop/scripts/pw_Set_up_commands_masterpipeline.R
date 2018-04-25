#Set up all the commands we need to do the allignment
#Set up a job array

#one command
#cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --output _trimmed_R1.fastq.gz --paired-output _trimmed_R2.fastq.gz  _R1.fastq.gz _R2.fastq.gz

#set working directory to source file location

library(stringr)

base.cmd <- "cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -b A{100} -b T{100} -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -B A{100} -B T{100} --minimum-length 50 --trim-n --output"

load(file="filelist.rda")

filelist <- gsub("_R1_001.fastq.gz", "", filelist, fixed=TRUE)
filelist <- gsub("_R2_001.fastq.gz", "", filelist, fixed=TRUE)
filelist <- unique(filelist)

#Trim
cat("", file="S01_Trim_LearnMemRNA.txt")
foldin <- "/group/kinglab/Patricka/Learn_Mem/rna-seq_data_raw/" 
foldout <- "/group/kinglab/Patricka/Learn_Mem/rna-seq_data_trimmed/" 

for(ff in filelist) 
  {
  cmd.cut <- paste(base.cmd," ",foldout,ff, "_trimmed_R1_001.fastq.gz --paired-output ",
                   foldout,ff, "_trimmed_R2_001.fastq.gz ",foldin,ff, "_R1_001.fastq.gz ",foldin,ff, "_R2_001.fastq.gz", sep="")
  cat(cmd.cut,"\n", file="S01_Trim_LearnMemRNA.txt",append=TRUE)
  
}

cat('', file="S02_QC_LearnMemRNA.txt")
foldin <- "/group/kinglab/Patricka/Learn_Mem/rna-seq_data_trimmed/" 
foldout <- "/group/kinglab/Patricka/Learn_Mem/fastqc_files/" 

base.cmd <- "fastqc -o" #../fastqc_files RAPiD-Genomics_HJYM3BBXX_MIZ_117501_P01_WB12_i5-509_i7-75_S577_L002_trimmed_R1_001.fastq.gz

for(ff in filelist) 
{
  cmd.cut1 <- paste(base.cmd," ",foldout," ",foldin,ff,"_trimmed_R1_001.fastq.gz", sep="")
  cmd.cut2 <- paste(base.cmd," ",foldout," ",foldin,ff,"_trimmed_R2_001.fastq.gz", sep="")
  
    cat(cmd.cut1,"\n", file="S02_QC_LearnMemRNA.txt",append=TRUE)
    cat(cmd.cut2,"\n", file="S02_QC_LearnMemRNA.txt",append=TRUE)
  
}


#Align
cat("", file="S03_Align_LearnMemRNA.txt")
foldin <- "/group/kinglab/Patricka/Learn_Mem/rna-seq_data_trimmed/" 
foldout <- "/group/kinglab/Patricka/Learn_Mem/rna-seq_data_aligned/" 

base.cmd <- "hisat2 --dta -x /group/kinglab/Patricka/base_pop/indexes/bdgp6_tran/genome_tran -1 /storage/hpc/group/kinglab/Patricka/Learn_Mem/rna-seq_data_trimmed/"

for(ff in filelist) 
{
  cmd.cut <- paste(base.cmd,ff,"trimmed_R1_001.fastq.gz -2 ",ff, "trimmed_R2_001.fastq.gz",
                              " -S ",foldout,ff,"_aligned_rna-seq.sam",sep="")
  cat(cmd.cut,"\n", file="S03_Align_LearnMemRNA.txt",append=TRUE)
}

#Sam to Bam
cat("", file="S04_SamToBam_LearnMemRNA.txt")
foldin <- "/group/kinglab/Patricka/Learn_Mem/rna-seq_data_aligned/" 
foldout <- "/group/kinglab/Patricka/Learn_Mem/rna-seq_data_bam/" 

base.cmd <- "samtools sort -@ 8 -o "

for(ff in filelist) 
{
  cmd.cut <- paste(base.cmd,foldout,ff,"_aligned_rna-seq.bam ",foldin,ff,"_aligned_rna-seq.sam",sep="")
  cat(cmd.cut,"\n", file="S04_SamToBam_LearnMemRNA.txt",append=TRUE)
}

#Merge Bam files from multiple lanes
cat("", file="S05_MergeBam_LearnMemRNA.txt")
foldin <- "/group/kinglab/Patricka/Learn_Mem/rna-seq_data_bam/" 
foldout <- "/group/kinglab/Patricka/Learn_Mem/rna-seq_data_bam/" 

ss <- read.csv(file="/home/pwilliams/DSPR/Scripts/rna-seq_processing_data/MIZ_117501_SampleSheet.csv",header=TRUE,stringsAsFactors = FALSE)
ss_fastq <- gsub("_R1_001.fastq", "", ss$Sequence_Name, fixed=TRUE)
issues<-which(!(ss_fastq %in% filelist))

ss<-ss[-issues,]

#samtools merge -r C-3_O_S9_all.bam C-3_O_S9_R1_001_run2.bam C-3_O_S9_R1_001_run3.bam
base.cmd <- "samtools merge -rf "

sample.ids <-unique(ss$Customer_Code)

for(ff in sample.ids) 
{
  lanes <- subset(ss, Customer_Code==ff)
  rr <- gsub("_R1_001.fastq", "", lanes$Sequence_Name, fixed=TRUE)
  #foldout,ff,"_aligned_rna-seq.bam
  
  cmd.cut <- paste(base.cmd,foldout,ff,"_merged.bam ",foldin,rr[1],"_aligned_rna-seq.bam ",foldin,rr[2],"_aligned_rna-seq.bam ",sep="")
  cat(cmd.cut,"\n", file="S05_MergeBam_LearnMemRNA.txt",append=TRUE)
}


#StringTie
cat("", file="S06_Assemble_LearnMemRNA.txt")
cat("", file="assembled_file_list.txt")

foldin <- "/group/kinglab/Patricka/Learn_Mem/rna-seq_data_bam/" 
foldout <- "/group/kinglab/Patricka/Learn_Mem/rna-seq_data_assemble/" 

#where is gtf?
#flybase
#uc santa cruz
base.cmd <- "stringtie -p 8 -G /group/kinglab/Patricka/base_pop/genes/dmel-all-r6.18.gtf -o "

for(ff in sample.ids) 
{
  cmd.cut <- paste(base.cmd,foldout,ff,"_assembled.gtf ",foldin,ff,"_merged.bam",sep="")
  cat(cmd.cut,"\n", file="S06_Assemble_LearnMemRNA.txt",append=TRUE)
  #merge  
  
  cat(paste(foldout,ff,"_assembled.gtf",sep=""), "\n",file="assembled_file_list.txt",append=TRUE)
}

#
#
#
##talk to libby about this step


#Examine how the transcripts compare to the reference annotation
#gffcompare –r chrX_data/genes/chrX.gtf –G –o merged stringtie_merged.gtf
#gffcompare



#cat("", file="S07_CompareRef_LearnMemRNA.txt")
#cat("", file="compare_reference_genome.txt")

#foldin <- "/group/kinglab/Patricka/Learn_Mem/rna-seq_LearnMem_merged/" 
#foldout <- "/group/kinglab/Patricka/Learn_Mem/rna-seq_data_compareRef/" 


#base.cmd <- "gffcompare –r /group/kinglab/Patricka/base_pop/genes/dmel-all-r6.18.gtf -o "

#for(ff in sample.ids) 
#{
  #cmd.cut <- paste(base.cmd,foldout,ff,"_assembled.gtf ",foldin,ff,"_merged.bam",sep="")
  #cat(cmd.cut,"\n", file="S06_Assemble_LearnMemRNA.txt",append=TRUE)
  #merge  
  
  #cat(paste(foldout,ff,"_assembled.gtf",sep=""), "\n",file="assembled_file_list.txt",append=TRUE)
#}


#Estimate transcript abundances and create table counts for Ballgown
# stringtie –e –B -p 8 -G stringtie_merged.gtf -o ballgown/ERR188044/ERR188044_chrX.gtf ERR188044_chrX.bam

#stringtie
cat("", file="S08_Abundances_LearnMemRNA.txt")


foldin <- "/group/kinglab/Patricka/Learn_Mem/rna-seq_data_bam/" 
foldout <- "/group/kinglab/Patricka/Learn_Mem/rna-seq_data_ballgown/" 


base.cmd <- "stringtie -e -B -p 8 -G /storage/hpc/group/kinglab/Patricka/Learn_Mem/rna-seq_data_bam/rna-seq_LearnMem_merged.gft -o "

for(ff in sample.ids) 
{
  cmd.cut <- paste(base.cmd,foldout,ff,"_ballgown.gtf ", foldin,ff,"_merged.bam",sep="")
  cat(cmd.cut,"\n", file="S08_Abundances_LearnMemRNA.txt",append=TRUE)
   
  
}

#make csv

#see above for where we make sample.ids

samp.desc <- data.frame('id'=character(length=length(sample.ids)),
                        'treatment'=character(length=length(sample.ids)),
                        'replicate'=character(length=length(sample.ids)),stringsAsFactors = FALSE)

samp.desc$id <- sample.ids 
samp.desc$treatment <- substr(samp.desc$id, 1,(nchar(samp.desc$id)-2))
samp.desc$replicate <- substr(samp.desc$id, (nchar(samp.desc$id)),(nchar(samp.desc$id)))

write.csv(samp.desc, file="sample_description.csv", row.names=FALSE, quote=FALSE)

#Ballgown step in the tuxedo pipeline,  
#load R version R version 3.2.2
#load R packagees:
library(ballgown)
library(genefilter)
library(dplyr)
library(devtools)

#Load the phenotype data(.cvs file)

pheno_data <- read.csv("sample_description.csv")


#Read in the expression data that was calculated by StringTie.

expression_data <- ballgown(dataDir = "ballgown", samplePattern = "ERR", pData=pheno_data)


#Filter to remove low-abundance genes.


#ballgown
#cat("", file="S09_Abundances_LearnMemRNA.txt")


#foldin <- "/group/kinglab/Patricka/Learn_Mem/rna-seq_data_bam/" 
#foldout <- "/group/kinglab/Patricka/Learn_Mem/rna-seq_data_ballgown/" 


#base.cmd <- "stringtie  –e –B -p 8 -G /storage/hpc/group/kinglab/Patricka/Learn_Mem/rna-seq_data_bam/rna-seq_LearnMem_merged.gft -o "

#for(ff in sample.ids) 
#{
  #cmd.cut <- paste(base.cmd,foldin,ff,"_ballgown.gtf ", foldout,ff,"_merged.bam",sep="")
 #cat(cmd.cut,"\n", file="S08_Abundances_LearnMemRNA.txt",append=TRUE)

#}



