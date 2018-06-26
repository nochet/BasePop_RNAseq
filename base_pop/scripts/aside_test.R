
# read table prepped from the ballgown `bg_ballG_all_results.Rda` object

basep_all <- read.csv("../processed/results/ballG_all_results/logGene_associated_texpr.csv",
                      stringsAsFactors = FALSE)

phenDat <- read.csv("../processed/describe_samples_batch.csv", stringsAsFactors = FALSE)

# Set up data for sva

# sample data 
phenDat = phenDat 

# Only expression variables needed in matrix form
basep_all = basep_all %>%
  select(-geneNames, -geneIDs) %>%
  select(tname, everything())

nams <- basep_all$tname
rownames(basep_all) = make.names(nams, unique=TRUE) 
basep_all$tname <- NULL

basep_all <- as.matrix(basep_all)

# Get a matrix of batch-corrected values from svaseq

# 1. Perform svaseq
group <- as.factor(rep(c("C", "DR", "HS"), each=18))

modd1 <- model.matrix(~group)
modd0 <- cbind(modd1[,1])

svaseq.dat <- svaseq(basep_all, modd1, modd0, n.sv=n.sv)

# 2. Extract corrected values from svaseq object
# a) can use function
# https://support.bioconductor.org/p/87508/
# https://support.bioconductor.org/p/47350/

cleaningY = function(y, mod, svaobj) {
  X = cbind(mod, svaobj$sv)
  Hat = solve(t(X)%*%X)%*%t(X)
  beta = (Hat%*%t(y))
  P = ncol(mod)
  cleany = y-t(as.matrix(X[,-c(1:P)])%*%beta[-c(1:P),]) 
  return(cleany)
}

svaseq.dat <- cleaningY(basep_all, mod = modd1, svaobj = svaseq.dat)
plot(svaseq.dat)

# b) can use sva & limma
svaseq_basep_all <- svaseq(basep_all, mod, mod0 , n.sv)
covv <- cbind(svaseq_basep_all$sv[,1:n.sv])
svaseq_basep_all2 <- removeBatchEffect(basep_all, covariates = covv)
plot(svaseq_basep_all2)






