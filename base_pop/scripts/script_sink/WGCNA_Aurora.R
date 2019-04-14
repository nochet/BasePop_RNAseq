# AURORA'S WGCNA CODE
# Preprint: Network analysis allows to unravel breast cancer molecular features and to identify novel targets

# Figure 1

#data are in matrix metabric, metadata are in data.frame metadata
#moduleColors is the vector with colors assigned to each gene by WGCNA

##calculate modules eigengenes
MEs_metabric= moduleEigengenes(t(metabric), moduleColors)$eigengenes

####selection of metadata columns
meta2<-c("grade", "age_at_diagnosis")

###module-trait correlations and p-values -- only 2 traits in this case
moduleTraitCor = list();
moduleTraitPvalue = list();
# Calculate the correlations
for (set in meta2){
  moduleTraitCor[[set]] = cor(MEs_metabric, as.numeric(as.character(meta[,set])), use = "p");
  moduleTraitPvalue[[set]] = corPvalueFisher(moduleTraitCor[[set]], ncol(metabric));
}

moduleTraitCor_tot<-moduleTraitCor[[1]]
moduleTraitPvalue_tot<-moduleTraitPvalue[[1]]

for (set in 2){ ##only 2 because there are only 2 traits, otherwise 2:n
  moduleTraitCor_tot<-cbind(moduleTraitCor_tot,moduleTraitCor[[set]])
  moduleTraitPvalue_tot<-cbind(moduleTraitPvalue_tot, moduleTraitPvalue[[set]])
}
colnames(moduleTraitCor_tot)<-meta2
colnames(moduleTraitPvalue_tot)<-meta2

#####discard "grey" module
MEs_metabric<-MEs_metabric[,-which(colnames(MEs_metabric)=="MEgrey")]

#### Convert numerical lables to colors for labeling of modules in the plot
MEColors = gsub("ME", "", colnames(MEs_metabric))
MEColorNames = alt_names[match(gsub( "ME", "",colnames(MEs_metabric)), labels2colors(1:20))];

#### Open a suitably sized window (the user should change the window size if necessary)
moduleTraitPvalue_tot[,1]<-p.adjust(moduleTraitPvalue_tot[,1])
moduleTraitPvalue_tot[,2]<-p.adjust(moduleTraitPvalue_tot[,2])

png("ModuleTraitRelationship_metabric.png", res = 450, width=2500, height = 4000)
# Plot the module-trait relationship table

textMatrix = paste(signif(moduleTraitCor_tot, 2), "\n(",
                   signif(moduleTraitPvalue_tot, 1), ")", sep = "");

###p-values less than 10^-16 are converted to 10^-16, the most significant p-value that can be reliably calculated
textMatrix[moduleTraitPvalue_tot<2.2*10^(-16)]<-paste(signif(moduleTraitCor_tot[moduleTraitPvalue_tot<2.2*10^(-16)], 2), "\n(",
                                                      "<2.2e-16)", sep = "");
dim(textMatrix) = dim(moduleTraitCor_tot)
par(mar = c(6, 14, 3, 2.2));
labeledHeatmap(Matrix = moduleTraitCor_tot,
               xLabels = meta2,
               yLabels = MEColorNames,
               ySymbols = MEColors,
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.8,
               cex.lab.x=1,
               cex.lab.y=1,
               zlim = c(-1,1),
               main = paste("Module - clinical feature"))
dev.off();







# Figure 5

####function to calculate the enrichment of each set of genes for module genes
#####genes -> set of genes to calculate the enrichment for
#####dataset -> microarray dataset used for the selection of differentially expressed genes (needed for the background)

fishertest_basal<-function(genes, dataset){
  fisherp<-matrix(ncol=length(unique(moduleColors_basal)), nrow=1)
  for(i in 1:length(unique(moduleColors_basal))){
    metabric_tmp<-metabric[rownames(metabric) %in% rownames(dataset),]
    moduleColors_basal_tmp<-moduleColors_basal[rownames(metabric) %in% rownames(dataset)]
    counts<-matrix(c(length(intersect(rownames(metabric_tmp)[which(moduleColors_basal_tmp %in% unique(moduleColors_basal)[i])], genes )),
                     length(intersect(rownames(metabric_tmp)[-which(moduleColors_basal_tmp %in% c("grey", unique(moduleColors_basal)[i]))], genes )),
                     length(which(moduleColors_basal_tmp %in% unique(moduleColors_basal)[i])),
                     nrow(metabric_tmp)-length(which(moduleColors_basal_tmp=="grey"))), nrow=2)
    fisherp[1,i]<-fisher.test(counts, alternative = "greater")[[1]]
  }
  
  colnames(fisherp)<-unique(moduleColors_basal)
  return(fisherp)
}


#####FOXM1
#import data
FOXM1_GSE2222=read.csv("GSE2222_series_matrix.txt", sep="\t", row.names = 1)

#gene id conversion
#the for cycle selects one probe per gene: if more probes map to the same gene, only the most expressed one is kept

library(hgu133a.db)
x <- hgu133aSYMBOL
mapped_probes <- mappedkeys(x)
xx <- as.list(x[mapped_probes])
vals <- sapply(xx, as.vector)
adf <- data.frame(probe=names(vals), gene=vals)

annotation_sel=adf[match( rownames(FOXM1_GSE2222), adf[,1]),2]

FOXM1_GSE2222<-FOXM1_GSE2222[which(is.na(annotation_sel)==F),]
annotation_sel<-na.omit(annotation_sel)
annotation_sel=as.character(annotation_sel)
values<-unique(annotation_sel)
for(i in 1:length(values)){
  if(length(which(annotation_sel==values[i]))>1){
    m=which.max(rowMeans(FOXM1_GSE2222[which(annotation_sel==values[i]),], na.rm=T))
    FOXM1_GSE2222=FOXM1_GSE2222[-which(annotation_sel==values[i])[-m],]
    annotation_sel=annotation_sel[-which(annotation_sel==values[i])[-m]]
  }
}
rownames(FOXM1_GSE2222)=annotation_sel

###t.test for the significant differential expression of each gene in the data matrix to retrieve downregulated genes
test_pval<-c()
for(i in 1:nrow(FOXM1_GSE2222)){
  if(length(unique(as.numeric(FOXM1_GSE2222[i,c(4:9)])))<3){
    test_pval<-c(test_pval, 1)
  } else{
    test_pval<-c(test_pval, t.test(FOXM1_GSE2222[i,4:6], FOXM1_GSE2222[i,7:9], alternative="greater")[[3]])
  }
}

###t.test for the significant differential expression of each gene in the data matrix to retrieve upregulated genes
test_pval_less<-c()
for(i in 1:nrow(FOXM1_GSE2222)){
  if(length(unique(as.numeric(FOXM1_GSE2222[i,c(4:9)])))<3){
    test_pval_less<-c(test_pval_less, 1)
  } else{
    test_pval_less<-c(test_pval_less, t.test(FOXM1_GSE2222[i,4:6], FOXM1_GSE2222[i,7:9], alternative="less")[[3]])
  }
}

#selection of DEGs with a pvalue<0.05
genes_up_FOXM1_GSE2222<-rownames(FOXM1_GSE2222)[which(test_pval_less<0.05)]
genes_down_FOXM1_GSE2222<-rownames(FOXM1_GSE2222)[which(test_pval<0.05)]

###calculate the enrichment for module genes
fishertest_basal(genes_up_FOXM1_GSE2222, FOXM1_GSE2222)
fishertest_basal(genes_down_FOXM1_GSE2222, FOXM1_GSE2222)


##########After having repeated the analysis for a series of datasets, the pvalues are plotted in a heatmap

DOWN<-fishertest_basal(genes_down_FOXM1_GSE2222, FOXM1_GSE2222)
DOWN<-rbind(DOWN, fishertest_basal(genes_down_FOXM1, FOXM1))
DOWN<-rbind(DOWN, fishertest_basal(genes_down_FOXM1_GSE25741, FOXM1_GSE25741))
DOWN<-rbind(DOWN, fishertest_basal(genes_down_PTTG1, PTTG1))
DOWN<-rbind(DOWN, fishertest_basal(genes_down_EZH2_GSE48979, EZH2_GSE48979)) 
DOWN<-rbind(DOWN, fishertest_basal(genes_down_EZH2_GSE36939_HCC70, EZH2_GSE36939))
DOWN<-rbind(DOWN, fishertest_basal(genes_down_EZH2_GSE36939_468, EZH2_GSE36939))
DOWN<-rbind(DOWN, fishertest_basal(genes_up_EZH2_GSE36939_OE, EZH2_GSE36939)) 
DOWN<-rbind(DOWN, fishertest_basal(genes_up_EZH2_GSE103242, EZH2_RPKM))
DOWN<-rbind(DOWN, fishertest_basal(genes_down_PARP1_GSE34817_10nM, PARP1_GSE34817)) 
DOWN<-rbind(DOWN, fishertest_basal(genes_down_PARP1_GSE34817_100nM, PARP1_GSE34817)) 

DOWN<-rbind(DOWN, fishertest_basal(genes_up_FOXC1, FOXC1_GSE73234))
DOWN<-rbind(DOWN, fishertest_basal(genes_down_HMGA1_GSE35525, HMGA1_GSE35525))
DOWN<-rbind(DOWN, fishertest_basal(genes_down_TCF7L1_GSE38893, TCF7L1_GSE38893)) 
DOWN<-rbind(DOWN, fishertest_basal(genes_down_SSRP1_T47D, SSRP1)) 
DOWN<-rbind(DOWN, fishertest_basal(genes_down_SSRP1_MCF7, SSRP1)) 
DOWN<-rbind(DOWN, fishertest_basal(genes_down_ELF5_GSE30405_HCC, ELF5_GSE30405))
DOWN<-rbind(DOWN, fishertest_basal(genes_down_ELF5_GSE30405_MCF7, ELF5_GSE30405)) 
DOWN<-rbind(DOWN, fishertest_basal(genes_down_ELF5_GSE30405_T47D, ELF5_GSE30405)) 


rownames(DOWN)<-c("FOXM1_GSE2222_BT-20", "FOXM1_GSE55204_MCF7", "FOXM1_GSE25741_MDA-MB-231","PTTG1_GSE48928_MCF7", "EZH2_GSE48979_MDA-MB-231","EZH2_GSE36939_HCC70", "EZH2_GSE36939_MDA468","EZH2_GSE36939_HCC70_OE", "EZH2_GSE103242_MCF7_OE",
                  "PARP1_GSE34817_10nM_MDA-MB-436","PARP1_GSE34817_100nM_MDA-MB-436",
                  "FOXC1_GSE73234_MDA-MB-231_OE", "HMGA1_GSE35525_MDA-MB-231", "TCF7L1_GSE38893_MDA-MB-468", "SSRP1_GSE92281_T47D", "SSRP1_GSE92281_MCF7", "ELF5_GSE30405_HCC", "ELF5_GSE30405_MCF7_OE","ELF5_GSE30405_T47D_OE" )

#assign alternative names to modules, names kept in the vector "alt_names_b"
colnames(DOWN)<-alt_names_b[match(colnames(DOWN), labels2colors(1:20))]
#remove "grey"
DOWN<-DOWN[,-1]
#prepare values to plot: -log10 of pvalues
logDOWN<-(-log10(DOWN))
#maximum significance reliable is 2.2*10-16
logDOWN[logDOWN>15.65758]<-15.65758

#two versions of the same heatmap scaling or not
library(pheatmap)

png("heatmap_TF_datasets_scaled.png",575,500)
pheatmap(logDOWN, scale="row")
dev.off()

png("heatmap_TF_datasets.png",575,500)
pheatmap(logDOWN)
dev.off()
