
library(DESeq2)

args <- commandArgs(TRUE)
Ni <- as.numeric(args[1])

RNs <- read.table(file="../../processed/RandNs.txt",header=TRUE)
cat(as.integer(RNs[Ni,1]))
set.seed(as.integer(RNs[Ni,1]))
 


## Simulation to see what fold changes we get by chance to look at similarity of diet responses

# Get obs data

# Expression matrix
countdata <- read.csv("../../processed/DESEQ/Expr_countData.csv")
rownames(countdata) <- countdata[,1]
countdata[,1] <- NULL
countdata <- as.matrix(countdata)

phenDat.sva <- read.csv("../../processed/DESEQ/sampleDat_with_SV.csv")
rownames(phenDat.sva) <- phenDat.sva[,1]
phenDat.sva[,1] <- NULL
phenDat.sva$batch <- as.factor(phenDat.sva$batch)
phenDat.sva$treat_tissue <- paste(phenDat.sva$treatment,"_",phenDat.sva$tissue, sep="")

phenDat.sva <- phenDat.sva[c(seq(1,54, by=3), seq(2,54, by=3), seq(3,54, by=3)),]
countdata <- countdata[,c(seq(1,54, by=3), seq(2,54, by=3), seq(3,54, by=3))]


#shuffle with 2 per treatment for each organ

Cs <- c(sample(seq(1,6)), sample(seq(19,24)), sample(seq(37,42)))
DRs <- c(sample(seq(7,12)), sample(seq(25,30)), sample(seq(43,48)))
HSs <- c(sample(seq(13,18)), sample(seq(31,36)), sample(seq(49,54)))

ss <- c(Cs[1:2],DRs[1:2],HSs[1:2],Cs[3:4],DRs[3:4],HSs[3:4],Cs[5:6],DRs[5:6],HSs[5:6],
        Cs[7:8],DRs[7:8],HSs[7:8],Cs[9:10],DRs[9:10],HSs[9:10],Cs[11:12],DRs[11:12],HSs[11:12],
        Cs[13:14],DRs[13:14],HSs[13:14],Cs[15:16],DRs[15:16],HSs[15:16],Cs[17:18],DRs[17:18],HSs[17:18]) 

#ss <- c(sample(seq(1,18)),sample(seq(19,36)),sample(seq(37,54)))

phenDat.sva.s <- phenDat.sva
phenDat.sva.s[,"batch"] <- phenDat.sva.s[ss,"batch"]
phenDat.sva.s[,"SV1"] <- phenDat.sva.s[ss,"SV1"]

countdata.s <- countdata[,ss]
colnames(countdata.s) <- colnames(countdata)

phenDat.sva[ss,]

dds <- DESeqDataSetFromMatrix(countData = countdata.s, 
                                   colData = phenDat.sva.s, 
                                   design = ~ treat_tissue)
dds <- estimateSizeFactors(dds)
nc <- counts(dds, normalized=TRUE)
filter <- rowSums(nc >= 10) >= 2
#dds <- dds[filter,]
dds <- estimateDispersions(dds, fitType = "local")
dds <- nbinomWaldTest(dds, maxit = 1000, useOptim = TRUE)


dds.01 <- DESeqDataSetFromMatrix(countData = countdata.s, 
                                   colData = phenDat.sva.s, 
                                   design = ~  tissue + treatment)
lrt.treatment <- DESeq(dds.01, test="LRT", reduced=~  tissue)
lrt.treatAll <- results(lrt.treatment)

all.equal(rownames(lrt.treatAll), rownames(dds))
all.equal(rownames(lrt.treatAll), names(filter))

lrt.treatAll <- lrt.treatAll[filter,]
dds <- dds[filter,]


organs <- c('B','H','O')
treats <- c('DR','C','HS')


all.fc.dat.s <- lrt.treatAll

for(oo in organs)
{

  FC_DRC <- lfcShrink(dds, contrast=c("treat_tissue",paste("DR_",oo,sep=""),paste("C_",oo,sep="")), type="normal")
  
  FC_HSC <- lfcShrink(dds, contrast=c("treat_tissue",paste("HS_",oo,sep=""),paste("C_",oo,sep="")), type="normal")
  
  FC_DRHS <- lfcShrink(dds, contrast=c("treat_tissue",paste("DR_",oo,sep=""),paste("HS_",oo,sep="")), type="normal")
  
  FC_CHS <- lfcShrink(dds, contrast=c("treat_tissue",paste("C_",oo,sep=""),paste("HS_",oo,sep="")), type="normal")
  
  FC_CDR <- lfcShrink(dds, contrast=c("treat_tissue",paste("C_",oo,sep=""),paste("DR_",oo,sep="")), type="normal")
  
  FC_HSDR <- lfcShrink(dds, contrast=c("treat_tissue",paste("HS_",oo,sep=""),paste("DR_",oo,sep="")), type="normal")
  
  
  
  colnames(FC_DRC) <- c('baseMean',
                        paste("FC_DR_",oo,sep=""),
                        paste("FCse_DR_",oo,sep=""),
                        paste("stat_DR_",oo,sep=""),
                        paste("p_DR_",oo,sep=""),
                        paste("padj_DR_",oo,sep=""))
  
  colnames(FC_HSC) <- c('baseMean',
                        paste("FC_HS_",oo,sep=""),
                        paste("FCse_HS_",oo,sep=""),
                        paste("stat_HS_",oo,sep=""),
                        paste("p_HS_",oo,sep=""),
                        paste("padj_HS_",oo,sep=""))
  
  colnames(FC_DRHS) <- c('baseMean',
                         paste("FC_DRHS_",oo,sep=""),
                         paste("FCse_DRHS_",oo,sep=""),
                         paste("stat_DRHS_",oo,sep=""),
                         paste("p_DRHS_",oo,sep=""),
                         paste("padj_DRHS_",oo,sep=""))
  
  colnames(FC_CHS) <- c('baseMean',
                        paste("FC_CHS_",oo,sep=""),
                        paste("FCse_CHS_",oo,sep=""),
                        paste("stat_CHS_",oo,sep=""),
                        paste("p_CHS_",oo,sep=""),
                        paste("padj_CHS_",oo,sep=""))
  
  colnames(FC_CDR) <- c('baseMean',
                        paste("FC_CDR_",oo,sep=""),
                        paste("FCse_CDR_",oo,sep=""),
                        paste("stat_CDR_",oo,sep=""),
                        paste("p_CDR_",oo,sep=""),
                        paste("padj_CDR_",oo,sep=""))
  
  colnames(FC_HSDR) <- c('baseMean',
                         paste("FC_HSDR_",oo,sep=""),
                         paste("FCse_HSDR_",oo,sep=""),
                         paste("stat_HSDR_",oo,sep=""),
                         paste("p_HSDR_",oo,sep=""),
                         paste("padj_HSDR_",oo,sep=""))
  
  
  all.fc.dat.s <- cbind(all.fc.dat.s, FC_DRC[,2:6],FC_HSC[,2:6],FC_DRHS[,2:6],
                      FC_CHS[,2:6],FC_CDR[,2:6],FC_HSDR[,2:6])
cat(oo,"\n")  
}

cor.dat<- c(cor(all.fc.dat.s[,"FC_DR_B"],all.fc.dat.s[,"FC_HS_B"]),
            cor(all.fc.dat.s[,"FC_HSDR_B"],all.fc.dat.s[,"FC_CDR_B"]),
            cor(all.fc.dat.s[,"FC_DRHS_B"],all.fc.dat.s[,"FC_CHS_B"]),
            cor(all.fc.dat.s[,"FC_DR_H"],all.fc.dat.s[,"FC_HS_H"]),
            cor(all.fc.dat.s[,"FC_HSDR_H"],all.fc.dat.s[,"FC_CDR_H"]),
            cor(all.fc.dat.s[,"FC_DRHS_H"],all.fc.dat.s[,"FC_CHS_H"]),
            cor(all.fc.dat.s[,"FC_DR_O"],all.fc.dat.s[,"FC_HS_O"]),
            cor(all.fc.dat.s[,"FC_HSDR_O"],all.fc.dat.s[,"FC_CDR_O"]),
            cor(all.fc.dat.s[,"FC_DRHS_O"],all.fc.dat.s[,"FC_CHS_O"]))

ns.dat<- c((length(which(all.fc.dat.s[,"FC_DR_B"] > 0 & all.fc.dat.s[,"FC_HS_B"] > 0)) + length(which(all.fc.dat.s[,"FC_DR_B"] < 0 & all.fc.dat.s[,"FC_HS_B"] < 0)))/nrow(all.fc.dat.s),
           (length(which(all.fc.dat.s[,"FC_HSDR_B"] > 0 & all.fc.dat.s[,"FC_CDR_B"] > 0)) + length(which(all.fc.dat.s[,"FC_HSDR_B"] < 0 & all.fc.dat.s[,"FC_CDR_B"] < 0)))/nrow(all.fc.dat.s),
           (length(which(all.fc.dat.s[,"FC_DRHS_B"] > 0 & all.fc.dat.s[,"FC_CHS_B"] > 0)) + length(which(all.fc.dat.s[,"FC_DRHS_B"] < 0 & all.fc.dat.s[,"FC_CHS_B"] < 0)))/nrow(all.fc.dat.s),
           (length(which(all.fc.dat.s[,"FC_DR_H"] > 0 & all.fc.dat.s[,"FC_HS_H"] > 0)) + length(which(all.fc.dat.s[,"FC_DR_H"] < 0 & all.fc.dat.s[,"FC_HS_H"] < 0)))/nrow(all.fc.dat.s),
           (length(which(all.fc.dat.s[,"FC_HSDR_H"] > 0 & all.fc.dat.s[,"FC_CDR_H"] > 0)) + length(which(all.fc.dat.s[,"FC_HSDR_H"] < 0 & all.fc.dat.s[,"FC_CDR_H"] < 0)))/nrow(all.fc.dat.s),
           (length(which(all.fc.dat.s[,"FC_DRHS_H"] > 0 & all.fc.dat.s[,"FC_CHS_H"] > 0)) + length(which(all.fc.dat.s[,"FC_DRHS_H"] < 0 & all.fc.dat.s[,"FC_CHS_H"] < 0)))/nrow(all.fc.dat.s),
           (length(which(all.fc.dat.s[,"FC_DR_O"] > 0 & all.fc.dat.s[,"FC_HS_O"] > 0)) + length(which(all.fc.dat.s[,"FC_DR_O"] < 0 & all.fc.dat.s[,"FC_HS_O"] < 0)))/nrow(all.fc.dat.s),
           (length(which(all.fc.dat.s[,"FC_HSDR_O"] > 0 & all.fc.dat.s[,"FC_CDR_O"] > 0)) + length(which(all.fc.dat.s[,"FC_HSDR_O"] < 0 & all.fc.dat.s[,"FC_CDR_O"] < 0)))/nrow(all.fc.dat.s),
           (length(which(all.fc.dat.s[,"FC_DRHS_O"] > 0 & all.fc.dat.s[,"FC_CHS_O"] > 0)) + length(which(all.fc.dat.s[,"FC_DRHS_O"] < 0 & all.fc.dat.s[,"FC_CHS_O"] < 0)))/nrow(all.fc.dat.s))

all.fc.dat.sig <- all.fc.dat.s[is.na(all.fc.dat.s$padj)==FALSE,]
all.fc.dat.sig <- all.fc.dat.s[all.fc.dat.sig$padj<=0.05,]

if(nrow(all.fc.dat.sig)>0){
cor.dat.sig<- c(cor(all.fc.dat.sig[,"FC_DR_B"],all.fc.dat.sig[,"FC_HS_B"]),
            cor(all.fc.dat.sig[,"FC_DR_H"],all.fc.dat.sig[,"FC_HS_H"]),
            cor(all.fc.dat.sig[,"FC_DR_O"],all.fc.dat.sig[,"FC_HS_O"]))

ns.dat.sig<- c((length(which(all.fc.dat.sig[,"FC_DR_B"] > 0 & all.fc.dat.sig[,"FC_HS_B"] > 0)) + length(which(all.fc.dat.sig[,"FC_DR_B"] < 0 & all.fc.dat.sig[,"FC_HS_B"] < 0)))/nrow(all.fc.dat.sig),
           (length(which(all.fc.dat.sig[,"FC_DR_H"] > 0 & all.fc.dat.sig[,"FC_HS_H"] > 0)) + length(which(all.fc.dat.sig[,"FC_DR_H"] < 0 & all.fc.dat.sig[,"FC_HS_H"] < 0)))/nrow(all.fc.dat.sig),
           (length(which(all.fc.dat.sig[,"FC_DR_O"] > 0 & all.fc.dat.sig[,"FC_HS_O"] > 0)) + length(which(all.fc.dat.sig[,"FC_DR_O"] < 0 & all.fc.dat.sig[,"FC_HS_O"] < 0)))/nrow(all.fc.dat.sig))
}else{
  cor.dat.sig <- NA
  ns.dat.sig <- NA
  
}

output <- list(cor.dat, ns.dat, cor.dat.sig, ns.dat.sig, nrow(all.fc.dat.sig))

output

save(output, file=paste("../../processed/FC_sim_",Ni,".rda", sep="" ))
