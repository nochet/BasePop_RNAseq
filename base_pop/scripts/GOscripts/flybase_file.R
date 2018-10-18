library(tidyverse)

#includes Dmel genes
GO_fly <- read_lines(file = "../../processed/gene_association.fb",
                                       skip = 5) %>% 
  str_split(pattern = "\t", simplify = TRUE) %>% 
  as_tibble()

colnames(GO_fly) <- c("DB","gname","g_symbol",
                      "Qualifier","GOid","reference",
                      "Evidence","sourceE","Aspect","full_gname",
                      "altname","type","taxon","date","sourceA")



