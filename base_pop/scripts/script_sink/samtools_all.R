#samtools sort -@ 8 -o C-3_O_S9_R1_001_run2.bam C-3_O_S9_R1_001_run2.sam
# samtools sort -o processed/HS-3_O_S45_R1_001_run3.bam processed/HS-3_O_S45_R1_001_run3.sam

ll<-list.files("../processed/", pattern="sam$")

for(ii in ll)
{
   cc<-strsplit(ii,".",fixed=TRUE)[[1]][1]
   cat(paste("samtools sort -@ 8 -o ../processed/",cc,".bam ", "../processed/", ii,sep=""), "\n", file=paste("samtools_all",".txt",sep=""),append=TRUE)

}

