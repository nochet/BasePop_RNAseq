library(DSPRqtl)
library(tidyverse)
data("positionlist_wgenetic")
load(file = "/home/kingeg/Projects/DSPRgeneral/Convert5_6/DSPR_5_6_parsed.rda")

source("mappingfunctions.R")

load(file="../Data/sig_ths.rda")#fdr.out
load(file="../Data/Lodscores_3traits.rda")#obs.LL

#th.2 <- fdr.out$fwer
th.l <- fdr.out$fdr$th[min(which(fdr.out$fdr$fdr<=0.05))]



#find peaks
peak.i <- apply(obs.LL[,1:3],2,function(x) peakInfo(x,cM=poslist[,c('chr','Gpos')] ,th=th.l, tol.dist=3))



ci.peak <- vector(mode="list", length=3)
names(ci.peak)<-names(peak.i)

for(kk in 1:3){
qq <- poslist[,1:3]
qq$LOD <- obs.LL[,kk]
pp.set <- peak.i[[kk]]
pp.set$lp <- NA
pp.set$up <- NA
pp.set$lg <- NA
pp.set$ug <- NA
pp.set$ulod <- NA
pp.set$llod <- NA
pp.set$chrR6 <- NA
pp.set$PposR6 <- NA
pp.set$lpR6 <- NA
pp.set$upR6 <- NA



for(jj in 1:nrow(pp.set))
{
ff<-findCI(pp.set$chr[jj],pp.set$Ppos[jj], qtldat=qq, method='BCI')
pp.set$up[jj] <- ff$Ppos[2]
pp.set$lp[jj] <- ff$Ppos[1]
pp.set$ug[jj] <- ff$Gpos[2]
pp.set$lg[jj] <- ff$Gpos[1]
pp.set$ulod[jj] <- ff$LOD[2]
pp.set$llod[jj] <- ff$LOD[1] 
pp.set$chrR6[jj] <- coord.table[which(coord.table$R5chr==pp.set$chr[jj] & coord.table$R5pos==pp.set$Ppos[jj]),'R6chr']
pp.set$PposR6[jj] <- as.numeric(coord.table[which(coord.table$R5chr==pp.set$chr[jj] & coord.table$R5pos==pp.set$Ppos[jj]),'R6pos'])

pp.set$lpR6[jj] <- as.numeric(coord.table[which(coord.table$R5chr==pp.set$chr[jj] & coord.table$R5pos==ff$Ppos[1]),'R6pos'])
pp.set$upR6[jj] <- as.numeric(coord.table[which(coord.table$R5chr==pp.set$chr[jj] & coord.table$R5pos==ff$Ppos[2]),'R6pos'])

} 



pp.set<-pp.set[order(pp.set$chr,pp.set$Ppos),]

ci.peak[[kk]]<-pp.set
}


#clean up redundant peaks by hand 
#EXAMINE ON PLOT
#[[1]] 3R 17640000 
#[[2]] 2R  3130000 
#[[2]] 3L 15650000
#[[2]] 3L 23440000
#[[2]] 3R   820000

#[[3]] 2R  2720000 BCI has problem
#[[3]] 2R  7150000 
#[[3]] 3L 16740000 BCI has problem
#[[3]] 3R 13190000 BCI has problem
#3R 15510000
#3R 17960000
#3R 20920000
# 3R 23450000
#re map with correction for large peak
#for now, consider these:
#2L 20110000
#3L  9380000
#3R 16810000



ci.peak[[1]] <- ci.peak[[1]][-which(ci.peak[[1]]$chr=='3R' & ci.peak[[1]]$Ppos==17640000),]

set.s <- c(which(ci.peak[[2]]$chr=='2R' & ci.peak[[2]]$Ppos==3130000),
           which(ci.peak[[2]]$chr=='3L' & ci.peak[[2]]$Ppos==15650000),
           which(ci.peak[[2]]$chr=='3L' & ci.peak[[2]]$Ppos==23440000),
           which(ci.peak[[2]]$chr=='3R' & ci.peak[[2]]$Ppos==820000))
ci.peak[[2]] <- ci.peak[[2]][-set.s,]

set.s <- c(which(ci.peak[[3]]$chr=='2L' & ci.peak[[3]]$Ppos==20110000),
           which(ci.peak[[3]]$chr=='3L' & ci.peak[[3]]$Ppos==9380000),
           which(ci.peak[[3]]$chr=='3R' & ci.peak[[3]]$Ppos==16810000))
ci.peak[[3]] <- ci.peak[[2]][set.s,]


save(ci.peak,file="../Data/Peaks_wCIs.rda")


###########


#load in all the datasets

#significant qtl peaks
load(file ="../Data/Peaks_wCIs.rda")
str(ci.peak)

#DE genes for Learning
load(file = "../Data/LearnresSVOrder.Rda")
str(LearnresSVorder)
LearnresSVorder$FBgn <- rownames(LearnresSVorder)
colnames(LearnresSVorder)

Learn_sig_genes <- as.data.frame(LearnresSVorder)

#colnames(LearnresSVorder)[which(names(LearnresSVorder) == "rownames")] <- "FBgn"

#DE genes for Memory
load(file = "../Data/MemresSVOrder.Rda")
str(MemresSVorder)
MemresSVorder$FBgn <- rownames(MemresSVorder)
colnames(MemresSVorder)

Mem_sig_genes <- as.data.frame(MemresSVorder)

#DE genes thermal Tolerance 

#interaction
load(file="../Data/TTlrt_inter.Rda")
str(TTlrt_inter)
TTlrt_inter$FBgn <- rownames(TTlrt_inter)
colnames(TTlrt_inter)

ThermTol_inter_sig_genes <- as.data.frame(TTlrt_inter)

#condition
load(file="../Data/TTlrt_condition.Rda")
str(TTlrt_condition)
TTlrt_condition$FBgn <- rownames(TTlrt_condition)
colnames(TTlrt_condition)

ThermTol_condition_sig_genes <- as.data.frame(TTlrt_condition)


#pool
load(file="../Data/TTlrt_pool.Rda")
str(TTlrt_pool)
TTlrt_pool$FBgn <- rownames(TTlrt_pool)
colnames(TTlrt_pool)

ThermTol_pool_sig_genes <- as.data.frame(TTlrt_pool)


#gene list from Fly Base

gene_map_table <- read_lines(file = "../Data/gene_map_table_fb_2018_04.tsv",
                skip = 6) %>% 
  str_split(pattern = "\t", simplify = TRUE) %>% 
  as_tibble() %>% 
  filter(V1 == "Dmel")

colnames(gene_map_table) <- c('spp','gname','FBgn','v3','cyt','pos')



gene_map_table$chr <- str_split(gene_map_table$pos, ":",simplify=TRUE)[,1]
temp.s1 <- str_split(gene_map_table$pos, ":",simplify=TRUE)[,2] %>% 
  str_split(fixed(".."),simplify=TRUE)
gene_map_table$startp <- as.numeric(temp.s1[,1])

gene_map_table$stopp <- as.numeric(str_split(temp.s1[,2],fixed("("),simplify = TRUE)[,1])

str(gene_map_table)
colnames(gene_map_table)


#gene_map_table <- as.numeric(gene_map_table$startp)

#merge DE genes with gene list from Fly base (merge by FBgn num)
Learn_gene_list_merged <- merge(Learn_sig_genes, gene_map_table, by="FBgn")
Mem_gene_list_merged <- merge(Mem_sig_genes, gene_map_table, by="FBgn")
ThermTol_inter_gene_list_merged <- merge(ThermTol_inter_sig_genes, gene_map_table, by="FBgn")
ThermTol_pool_gene_list_merged <- merge(ThermTol_pool_sig_genes, gene_map_table, by="FBgn")
ThermTol_condition_gene_list_merged <- merge(ThermTol_condition_sig_genes, gene_map_table, by="FBgn")



str(Learn_gene_list_merged)
colnames(Learn_gene_list_merged)


#select specific columns

Learn_gene_list_sub <- Learn_gene_list_merged[c(1,3,6,7,13:15)]
colnames(Learn_gene_list_sub)
Learn_gene_list_sub <- subset(Learn_gene_list_sub, padj <=0.05)

Mem_gene_list_sub <- Mem_gene_list_merged[c(1,3,6,7,13:15)]
colnames(Mem_gene_list_sub)
Mem_gene_list_sub <- subset(Mem_gene_list_sub, padj <=0.05)

ThermTol_inter_gene_list_sub <- ThermTol_inter_gene_list_merged[c(1,3,6,7,13:15)]
colnames(ThermTol_inter_gene_list_sub)
ThermTol_inter_gene_list_sub <- subset(ThermTol_inter_gene_list_sub, padj <=0.05)

ThermTol_pool_gene_list_sub <- ThermTol_pool_gene_list_merged[c(1,3,6,7,13:15)]
colnames(ThermTol_pool_gene_list_sub)
ThermTol_pool_gene_list_sub <- subset(ThermTol_pool_gene_list_sub, padj <=0.05)

ThermTol_condition_gene_list_sub <- ThermTol_condition_gene_list_merged[c(1,3,6,7,13:15)]
colnames(ThermTol_condition_gene_list_sub)
ThermTol_conditiion_gene_list_sub <- subset(ThermTol_condition_gene_list_sub, padj <=0.05)


#filter dataset to narrow down genes under the qtl peaks
#USE R6

#Learning
foc.peak <- ci.peak[[1]][1,]

learn.list <-vector(mode='list', length=nrow(ci.peak$LearnPowerTrans))

Learn_totdf<-data.frame(NULL)
for(row in 1:nrow(ci.peak$LearnPowerTrans)) 
  
  {
  
  foc.peak<-ci.peak[['LearnPowerTrans']][row,]
  lw <- (which(Learn_gene_list_sub$chr==foc.peak$chr & 
                 ((Learn_gene_list_sub$startp <= foc.peak$upR6 & Learn_gene_list_sub$stopp >= foc.peak$lpR6))))  
  
  Learn_genes_under_peak<-Learn_gene_list_sub[lw,]
  learn.list[[row]] <-Learn_genes_under_peak
  
 }

save(Learn_genes_under_peak,file="../Data/Learn_genes_under_peak.rda")


#Memory 

Mem_foc_peak <- ci.peak[["Memory_Mean"]][1,]

Mem.list <-vector(mode='list', length=nrow(ci.peak$Memory_Mean))

Mem_totdf<-data.frame(NULL)
for(row in 1:nrow(ci.peak$Memory_Mean)) 
  
  {
  
  foc.peak<-ci.peak[['Memory_Mean']][row,]
  mw <- (which(Mem_gene_list_sub$chr==foc.peak$chr & 
                 ((Mem_gene_list_sub$startp <= foc.peak$upR6 & Mem_gene_list_sub$stopp >= foc.peak$lpR6))))  
  
 
  Mem_Mem_genes_under_peak<-Mem_gene_list_sub[mw,]
  
  Mem.list[[row]] <-Mem_genes_under_peak
  
}


Mem_genes_under_peak

save(Mem_genes_under_peak,file="../Data/Mem_genes_under_peak.rda")


#Learning (for a single peak)
foc.peak <- ci.peak[['LearnPowerTrans']][1,]

lw <- (which(Learn_gene_list_sub$chr==foc.peak$chr & 
               ((Learn_gene_list_sub$startp <= foc.peak$upR6 & Learn_gene_list_sub$stopp >= foc.peak$lpR6))))  

Learn_gene_list_sub[lw,]




#make a plot showing locations & pvalues?

