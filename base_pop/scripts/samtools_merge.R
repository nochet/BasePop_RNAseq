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

