#load libraries. You might  need to install these. Remember install.packages().
library(tidyverse)
library(cowplot)
library(wesanderson)

#load the data files. Remember to change the file paths if you need to or set your working directory.
load(file="results_genes.rda")
load(file="mm.rda")
load(file="glfpkm.rda")
load(file="tlfpkm.rda")
load(file='gid.rda')

#Get transcript abundances for boxplot
tt<- as.data.frame(tlfpkm) %>% gather(id, fpkm)
tt$treatment <- rep(c("C","DR","HS"), each=4*nrow(tlfpkm))

p1 <- ggplot(tt, aes(x=id, y=fpkm,fill=treatment)) +
  geom_boxplot() +
  ylab(expression("log"[2]*"(FPKM + 1)"))
p1
#ggsave(p1, filename="FPKM_overall.pdf", width=6, height=4)
#the above function, ggsave will save this as a pdf to whereever you want. 

#get fold changes
foldc <- mm %>% gather(treat, foldchange, -id, -pval)
foldc <- foldc[-which(foldc$foldchange>10),]

foldc$sig <- rep('a',nrow(foldc))

#pick significance level. this is our level after correcting for multiple tests. 
#we can be more liberal if you want for the plot
foldc$sig[foldc$pval<=0.00006]<-'b'

#split up the signficant and not for plotting
foldc1<-subset(foldc, sig=='a')
foldc2<-subset(foldc, sig=='b')

#get colors for significant points from wes anderson color palettes. 
#if you wanted to set your own set of 10, just make a vector of a list of 10 colors.
# if we are more liberal to show more genes, we need more colors
ccs<-c(wes_palette('Royal1'), wes_palette('Moonrise3'),wes_palette("Cavalcanti"))[1:10]

p2 <- ggplot(foldc1, aes(x=treat, y=foldchange)) +
  geom_point(pch=16, position=position_jitter(width = 0.3),alpha=0.2,color='grey') +
  geom_boxplot(outlier.shape = NA, alpha=0.8, color='grey40') +
  geom_point(data=foldc2, aes(x=treat, y=foldchange), color=rep(ccs,3),position=position_jitter(width = 0.3))+
  ylab("Fold Change") +
  xlab("Treatments Compared")
p2

all.equal(foldc2$id[1:10], as.factor(gid$id))

legd<-data.frame('y'=seq(1:10), 'x'=1,'cc'=ccs, 'id'=gid$gname)
p2L <- ggplot(legd, aes(x=x, y=y)) +
  geom_point(size=3, color=ccs) +
  geom_text(aes(label=as.character(id)),hjust=-0.2,vjust=0,fontface="italic") +
  labs(x="", y="") +
  xlim(c(0.95,1.2)) +
  theme_void()
pcomb <-plot_grid(p2,p2L, align = 'h',rel_widths = c(3, 1.4))
ggsave(pcomb, filename="../Plots/foldchange_overall.pdf", width=4, height=5)

#show significant results
ww <- which(results_genes$qval<=0.05)
rowMeans(glfpkm[ww,])
results_genes[ww,]

ff.sig <- merge(foldc, gid, by='id')

ff.sig

#make plot for one gene. 
ii<-1

for(ii in 1:10)
  {
pp<-data.frame('FPKM'=glfpkm[ww[ii],], 'Treatment' = rep(c('C','DR','HS'),each=4))

p3 <- ggplot(pp, aes(x=Treatment, y=FPKM)) +
  geom_point(position=position_jitter(width = 0.05),color='steelblue', alpha=0.7,size=3) +
  ylim(c(min(pp$FPKM)-1,max(pp$FPKM)+1)) +
  geom_text(label=gid$gname[ii], x=1, y=max(pp$FPKM)+1, fontface="italic")+
  ylab(expression("log"[2]*"(FPKM + 1)")) 
fname <- paste("../Plots/Gene_",gid$gname[ii], ".pdf", sep="")
ggsave(p3, filename=fname, width=4, height=4)
#p3
}

