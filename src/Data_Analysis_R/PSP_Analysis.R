library(tidyverse) 
library(data.table)

modified_residue <- function(mod_rsd){
  residue = substr(mod_rsd,1,1)
  return(residue)
}

psp_db <- data.table::fread("/Users/aparnasrinivasan/Desktop/MScThesis/Srinivasan_MS_Thesis/Data/Raw_Data/PhosphoSitePlus_Dataset/Phosphorylation_site_dataset.tsv") %>%
  filter(ORGANISM == "human") %>% 
  mutate(residue = modified_residue(MOD_RSD))

psp_db %>% group_by(PROTEIN, MW_kD) %>% 
  #filter(MS_LIT > 5) %>% 
  summarise(n_sites = n_distinct(MOD_RSD)) %>% ungroup() -> db

summary(db$MW_kD) 

iqr = IQR(db$MW_kD)
third_quartile = quantile(db$MW_kD, .75)
mw_range = third_quartile  + iqr 

print(mw_range)

# Outliers removed 

ggplot(db %>% filter(MW_kD < mw_range), aes(x=n_sites)) + geom_histogram(binwidth = 1) + 
  xlab("Unique phosphorylation sites reported") + 
  ylab("Proteins") + 
  geom_vline(aes(xintercept=mean(n_sites)), linetype="dashed", color="red")


p1 <- ggplot(db, aes(x=n_sites)) + geom_histogram() + 
  xlab("Unique phosphorylation sites reported") + 
  ylab("Proteins") 

p2 <- ggplot(db %>% filter(n_sites < 50), aes(x=n_sites)) + geom_histogram(binwidth = 1) + 
 labs(x=NULL, y=NULL) +
  geom_vline(aes(xintercept=12.0965), linetype="dashed", color="red")


p1 + annotation_custom(ggplotGrob(p2), xmin= 125, xmax=650, ymin = 2500, ymax = 15000)



knitr::kable(db %>% summarise(mean(n_sites), median(n_sites)))

summary((db %>% filter(n_sites >12))$MW_kD)
summary((db %>% filter(n_sites >6))$MW_kD)

ggplot(psp_db, aes(x=residue)) + geom_bar() +
  xlab("Amino Acid") + 
  ylab("Unique phosphorylation sites")

