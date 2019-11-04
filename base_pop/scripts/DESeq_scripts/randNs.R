set.seed(39426124)

tt <- as.integer(runif(100,1,1000000))

write.table(tt, file="../../processed/RandNs.txt",sep="\t",row.names = FALSE)


for(i in 1:100)
{
pp <- paste("Rscript FoldChange_sim.R ",i," >temp",i,".txt 2>error",i,".txt &",sep="")  
cat(pp, "\n", file="../../processed/sims_cmds.txt",append=TRUE)  
}



