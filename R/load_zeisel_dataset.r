load_zeisel_dataset <- function(){
  # Data from: http://mousebrain.org/downloads.html
  library(tidyverse)
  library("rhdf5")
  library("snow")

  linURL = "https://storage.googleapis.com/linnarsson-lab-loom/l5_all.agg.loom"
  path = "l5_all.agg.loom"
  download.file(linURL, destfile=path)
  
  h5f <- H5Fopen(path)
  exp <- as.data.frame(t(h5f$matrix))
  exp$Gene <- h5f$row_attrs$Gene

  # Only keep genes with a unique name and tidy data
  exp <- exp %>% add_count(Gene) %>%
    filter(n==1) %>%
    select(-n) %>%
    gather(key = column,value=Expr,-Gene) %>%
    as.tibble()

  m2h <- read_tsv("Data/m2h.txt",col_types = "iccccc") %>%
    select(musName,entrez) %>%
    rename(Gene=musName) %>% rename(ENTREZ=entrez)

  cell_types <- cbind(column=as.character(paste0("V",1:265)),
                      Lvl1=h5f$col_attrs$TaxonomyRank1,
                      Lvl2=h5f$col_attrs$TaxonomyRank2,
                      Lvl3=h5f$col_attrs$TaxonomyRank3,
                      Lvl4=h5f$col_attrs$TaxonomyRank4,
                      Lvl5=h5f$col_attrs$ClusterName,
                      Description=h5f$col_attrs$Description,
                      NCells=h5f$col_attrs$NCells) %>%
    as.tibble() %>%
    mutate(NCells=as.numeric(NCells))

  exp_lvl5 <- inner_join(exp,cell_types,by="column") %>% ungroup() %>% rename(Expr_sum_mean=Expr)

  # Get summary stats

  sumstats_lvl5 <- exp_lvl5 %>%
    group_by(Lvl5) %>%
    summarise(mean_UMI=sum(Expr_sum_mean),
              Ngenes=sum(Expr_sum_mean>0))
  cell_types_arranged <- cell_types %>%
    group_by(Lvl5) %>%
    summarise(NCells=sum(NCells))
  sumstats_lvl5 <- inner_join(sumstats_lvl5,cell_types_arranged) %>%
    mutate(total_UMI=mean_UMI*NCells)

  # Cell types should have at least 200k reads.
  filter_lvl5 <- sumstats_lvl5 %>% filter(total_UMI>200000)
  exp_lvl5 <- filter(exp_lvl5,Lvl5%in%filter_lvl5$Lvl5)

  # Remove not expressed genes
  not_expressed <- exp_lvl5 %>%
    group_by(Gene) %>%
    summarise(total_sum=sum(Expr_sum_mean)) %>%
    filter(total_sum==0) %>%
    select(Gene) %>% unique()
  exp_lvl5 <- filter(exp_lvl5,!Gene%in%not_expressed$Gene)

  # Scale expression
  exp_lvl5 <- exp_lvl5 %>%
    group_by(Lvl5) %>%
    mutate(Expr_sum_mean=Expr_sum_mean*1e6/sum(Expr_sum_mean))

  exp_lvl5_wtSpec <- exp_lvl5 %>%
    group_by(Gene) %>%
    mutate(specificity=Expr_sum_mean/sum(Expr_sum_mean),
           max=max(Expr_sum_mean))

  file.remove(path)
  return(exp_lvl5_wtSpec)
}
