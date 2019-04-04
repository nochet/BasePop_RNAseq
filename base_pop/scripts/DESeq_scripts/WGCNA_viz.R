
phenDat.sva <- read.csv("../../processed/DESEQ/sampleDat_with_SV.csv")
eg.genes <- read.table(file="../../processed/DESEQ/Coexpression/WCGNA_eigengenes.txt", sep="\t", header=TRUE, stringsAsFactors = FALSE) 

all.dat <- cbind(phenDat.sva, eg.genes)

ccs <- colnames(all.dat[,6:26])
all.plots <- vector(mode="list", length=length(ccs))
names(all.plots) <- ccs

for(cc in ccs)
{

pp<- ggplot(all.dat, aes_string(x="treatment", y=cc, color="tissue")) +
  geom_point(position = position_jitter(width = 0.15)) +
  stat_summary(fun.y = mean, geom = "point", size = 3, pch=2)

all.plots[[cc]] <- pp

}

all.plots[[5]]
