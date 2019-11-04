
library(dplyr)

bvl <- read.table("~/MyGithub/BasePop_RNAseq/base_pop/processed/align_rate.txt",
                  sep = "\t", stringsAsFactors = FALSE)
bl <- bvl %>% 
  mutate(Percent = str_split(V1, ":", simplify = TRUE)[, 2],
         scope = str_split(V1, " ", simplify = TRUE)[, 2])

# use sub
