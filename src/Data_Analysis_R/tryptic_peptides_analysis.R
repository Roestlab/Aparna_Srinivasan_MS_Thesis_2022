library(data.table)
library(tidyverse)

tryptic_peptides <- data.table::fread("/Users/aparnasrinivasan/Desktop/MScThesis/Srinivasan_MS_Thesis/Data/Data_Analysis_Results/Human_tryptic_peptides_analysis/human_phosphoproteome_build.tsv") %>% 
                    mutate(n_sty = str_count(peptide_sequence, "[STY]"), peptide_length = str_length(peptide_sequence))


tryptic_peptides %>% summarise(n_distinct(peptide_sequence))

mean(tryptic_peptides$n_sty)

ggplot(tryptic_peptides, aes(x=n_sty)) + geom_histogram(binwidth = 1, center=1) + 
  xlab("Phosphorylatable residues (STY)") + 
  ylab("Peptide sequences") + 
  ggtitle("Phosphorylatable residues on human tryptic peptide database") + 
  geom_vline(aes(xintercept=median(n_sty)), linetype="dashed")

