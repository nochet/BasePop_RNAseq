#! /bin/env R

# Step 6: StringTie - transcript abundances and table of counts

```{r}
cat("", file="S06_abundances_pl.txt")


indir <- "/group/kinglab/enoch/MyGithub/BasePop_RNAseq/base_pop/processed/dotbams/" 
outdir <- "/group/kinglab/enoch/MyGithub/BasePop_RNAseq/base_pop/scripts/assembly/" 


base.cmd <- "stringtie -e -B -p 8 -G /group/kinglab/enoch/MyGithub/BasePop_RNAseq/base_pop/processed/perl_postPrep_S04/stringtie_merged_pl.gtf -o "

base.cmd <- "stringtie -e -B -p 8 -G /group/kinglab/enoch/MyGithub/BasePop_RNAseq/base_pop/processed/perl_postPrep/stringtie_merged_pl.gtf -o "

for(ff in un.set) 
{
  cmd.cut <- paste(base.cmd,outdir,ff,"/",ff,"_ballgown.gtf ", indir,ff,"_merged.bam",sep="")
  cat(cmd.cut,"\n", file="S06_abundances_pl.txt",append=TRUE)
   
}
```
