
library(dplyr)

bvl <- read.table("~/MyGithub/BasePop_RNAseq/base_pop/processed/align_rate.txt",
                  sep = "\t", stringsAsFactors = FALSE)
bl <- bvl %>% 
  mutate(Percent = str_split(V1, ":", simplify = TRUE)[, 2]) %>%
  select(Percent) %>%
  mutate(scope = str_split(Percent, "%", simplify = TRUE)[, 1]) %>%
  select(scope) 

bl$scope <- as.numeric(bl$scope) 
bl <- filter(bl, scope > 0)

max(bl)
min(bl)
median(bl$scope)


# Function to calculate mode
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

getmode(bl$scope)

