### R code from vignette source 'rnaseqAnalysis.Rnw'

###################################################
### code chunk number 1: load data (eval = FALSE)
###################################################
## library(ReportingTools)
## data(mockRnaSeqData)


###################################################
### code chunk number 2: run_edgeR (eval = FALSE)
###################################################
## library(edgeR)
## conditions <- c(rep("case",3), rep("control", 3))
## d <- DGEList(counts = mockRnaSeqData, group = conditions)
## d <- calcNormFactors(d)
## d <- estimateCommonDisp(d)
## ## Get an edgeR object
## edgeR.de <- exactTest(d)


###################################################
### code chunk number 3: edgeR_report (eval = FALSE)
###################################################
## library(lattice)
## rep.theme <- reporting.theme()
## ## Change symbol colors in plots
## rep.theme$superpose.symbol$col <- c("blue", "red")
## rep.theme$superpose.symbol$fill <- c("blue", "red")
## lattice.options(default.theme = rep.theme)
## ## Publish a report of the top 10 genes with p-values < 0.05 and log-fold change > 2
## ## In this case, the plots contain the counts from mockRnaSeqData, which are not normalized.
## ## The publish function does not normalize counts for the countTable argument to allow for
## ## flexibility in plotting various units (e.g. RPKM instead of counts).
## 
## deReport <- HTMLReport(shortName = 'RNAseq_analysis_with_edgeR',
##     title = 'RNA-seq analysis of differential expression using edgeR',
##     reportDirectory = "./reports")
## publish(edgeR.de, deReport, countTable=mockRnaSeqData,
## 	conditions=conditions, annotation.db = 'org.Mm.eg', 
## 	pvalueCutoff = .05, lfc = 2, n = 10, name="edgeR")
## finish(deReport)
## 
## ## If you would like to plot normalized counts, run the following commands instead:
## ## mockRnaSeqData.norm <- d$pseudo.counts
## ## publish(edgeR.de, deReport, mockRnaSeqData.norm, 
## ##        conditions, annotation.db = 'org.Mm.eg', 
## ## 	  pvalueCutoff = .05, lfc = 2, n = 10)
## ## finish(deReport)


###################################################
### code chunk number 4: edgeR_report (eval = FALSE)
###################################################
## d <- DGEList(counts = mockRnaSeqData, group = conditions)
## d <- calcNormFactors(d)
## design <- model.matrix(~conditions)
## d <- estimateGLMCommonDisp(d, design)
## d <- estimateGLMTrendedDisp(d, design)
## d <- estimateGLMTagwiseDisp(d, design)
## fit <- glmFit(d,design)
## edgeR.lrt <- glmLRT(fit, coef=2)
## 
## deReport2 <- HTMLReport(shortName = 'RNAseq_analysis_with_edgeR_2',
##     title = 'RNA-seq analysis of differential expression using edgeR (LRT)',
##     reportDirectory = "./reports")
## publish(edgeR.lrt, deReport2, countTable=mockRnaSeqData,
## 	conditions=conditions, annotation.db = 'org.Mm.eg', 
## 	pvalueCutoff = .05, lfc = 2, n = 10, name="edgeRlrt")
## finish(deReport2)


###################################################
### code chunk number 5: run_DESeq (eval = FALSE)
###################################################
## library(DESeq)
## cds<-newCountDataSet(mockRnaSeqData, conditions)
## cds<-estimateSizeFactors(cds)
## cds<-estimateDispersions(cds)
## res<-nbinomTest(cds,"control", "case" )


###################################################
### code chunk number 6: DESeq_report (eval = FALSE)
###################################################
## desReport <- HTMLReport(shortName = 'RNAseq_analysis_with_DESeq',
##     title = 'RNA-seq analysis of differential expression using DESeq',
##     reportDirectory = "./reports")
## publish(res,desReport,name="df",countTable=mockRnaSeqData, pvalueCutoff=0.05,
##     conditions=conditions,annotation.db="org.Mm.eg.db", 
##     expName="deseq",reportDir="./reports", .modifyDF=makeDESeqDF)
## finish(desReport)


###################################################
### code chunk number 7: run_DESeq2 (eval = FALSE)
###################################################
## library(DESeq2)
## conditions <- c(rep("case",3), rep("control", 3))
## mockRna.dse <- DESeqDataSetFromMatrix(countData = mockRnaSeqData,
##                         colData = as.data.frame(conditions), design = ~ conditions)
## colData(mockRna.dse)$conditions <- factor(colData(mockRna.dse)$conditions, levels=c("control", "case"))
## ## Get a DESeqDataSet object
## mockRna.dse <- DESeq(mockRna.dse)


###################################################
### code chunk number 8: DESeq2_report (eval = FALSE)
###################################################
## des2Report <- HTMLReport(shortName = 'RNAseq_analysis_with_DESeq2',
##     title = 'RNA-seq analysis of differential expression using DESeq2',
##     reportDirectory = "./reports")
## publish(mockRna.dse,des2Report, pvalueCutoff=0.05,
##     annotation.db="org.Mm.eg.db", factor = colData(mockRna.dse)$conditions,
##     reportDir="./reports")
## finish(des2Report)


###################################################
### code chunk number 9: Do GO analysis (eval = FALSE)
###################################################
## library(GOstats)
## library(org.Mm.eg.db)
## tt <- topTags(edgeR.de, n = 1000, adjust.method = 'BH', sort.by = 'p.value')
## selectedIDs <- rownames(tt$table)
## universeIDs <- rownames(mockRnaSeqData)
## goParams <- new("GOHyperGParams", 
##     geneIds = selectedIDs, 
##     universeGeneIds = universeIDs, 
##     annotation ="org.Mm.eg" , 
##     ontology = "MF", 
##     pvalueCutoff = 0.01,
##     conditional = TRUE, 
##     testDirection = "over")
## goResults <- hyperGTest(goParams)


###################################################
### code chunk number 10: make the GO report (eval = FALSE)
###################################################
## goReport <- HTMLReport(shortName = 'go_analysis_rnaseq',
## 	title = "GO analysis of mockRnaSeqData",
## 	reportDirectory = "./reports")
## publish(goResults, goReport, selectedIDs=selectedIDs, annotation.db="org.Mm.eg", 
## 	pvalueCutoff= 0.05)
## finish(goReport)


###################################################
### code chunk number 11: Do PFAM analysis (eval = FALSE)
###################################################
## library(Category)
## params <- new("PFAMHyperGParams", 
## 	geneIds= selectedIDs, 
## 	universeGeneIds=universeIDs, 
## 	annotation="org.Mm.eg",
## 	pvalueCutoff= 0.01,
## 	testDirection="over")
## PFAMResults <- hyperGTest(params)


###################################################
### code chunk number 12: make the PFAM report (eval = FALSE)
###################################################
## PFAMReport <- HTMLReport(shortName = 'pfam_analysis_rnaseq',
## 	title = "PFAM analysis of mockRnaSeqData",
## 	reportDirectory = "./reports")
## publish(PFAMResults, PFAMReport, selectedIDs=selectedIDs, annotation.db="org.Mm.eg",categorySize=5)
## finish(PFAMReport)


###################################################
### code chunk number 13: make the index page (eval = FALSE)
###################################################
## indexPage <- HTMLReport(shortName = "indexRNASeq",
##     title = "Analysis of mockRnaSeqData",
##     reportDirectory = "./reports")
## publish(Link(list(deReport,des2Report, goReport, PFAMReport), report = indexPage),
##     indexPage)
## finish(indexPage)


