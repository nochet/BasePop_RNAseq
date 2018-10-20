DEunderPeak <- function(DEgenes, arm, upper, lower)
{
  #DEgenes is a data frame with differentially expressed genes including chromosome (chr)
  # gene start (startp) and gene stop (stopp)
  
  #upper is the upper bound of the QTL CI in release 6 coordinates
  #lower is the lower bound of the QTL CI in release 6 coordinates
  
  lw <- (which(DEgenes$chr==arm & 
                 ((DEgenes$startp <= upper & DEgenes$stopp >= lower))))  
  
  return(DEgenes[lw,])
}



gene_map_table <- read_lines(file = "../../processed/gene_map_table_fb_2018_04.tsv",
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

write.csv(gene_map_table, "../../processed/DESEQ/gene_map_table.csv")

#merge with DEseq gene list
#subset to significant genes
#use DEunderpeak function


