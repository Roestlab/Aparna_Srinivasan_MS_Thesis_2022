library(data.table)
library(dplyr)
library(readr)
library(tidyverse)

ground_truth <- data.table::fread("/Users/aparnasrinivasan/Desktop/MScThesis/Srinivasan_MS_Thesis/Data/ground_truth_table.tsv")

diapasef_3_results <- data.table::fread("/Users/aparnasrinivasan/Desktop/MScThesis/Srinivasan_MS_Thesis/Data/Raw_Data_Search_Results/Pyprophet_Scored_Results/scored_results_dfs/diapasef_results_3_1.tsv") %>%
                      left_join(ground_truth %>% select(M_pool, Sequence) %>% rename(pool = M_pool) %>% mutate(mpool_gt = TRUE) %>% distinct(), by=c("pool", "Sequence"))

diapasef_3_results$mpool_gt = diapasef_3_results$mpool_gt %>% replace_na(FALSE)


library <- data.table::fread("/Users/aparnasrinivasan/Desktop/MScThesis/Srinivasan_MS_Thesis/Data/Data_Analysis_Results/Upool_spectral_library_generation/upool_library_formatted.tsv") %>% select(-pool)%>% rename(RT = iRT) %>% mutate(method= "Spectral library") %>% 
  left_join(ground_truth %>% select(M_pool, Sequence) %>% rename(pool = M_pool) %>% mutate(mpool_gt = TRUE) %>% distinct(), by=c("Sequence")) 

library$mpool_gt <- library$mpool_gt %>% replace_na(FALSE)

iRT_peptides <- unique((library %>% filter(mpool_gt == FALSE))$Sequence)



filtered_diapasef_3_results <- diapasef_3_results %>% 
                                filter(mpool_gt == TRUE) %>% 
                                filter(ms2_pep < 0.05) %>% 
                                filter(!grepl("v2", filename)) %>%
  left_join(library %>% select(FullPeptideName, Charge, IM) %>% distinct() %>% rename(library_IM = IM), by=c("FullPeptideName", "Charge"))# no data from the samples run in March 2022 for now 

write_tsv(filtered_diapasef_3_results, "/Users/aparnasrinivasan/Desktop/MScThesis/Srinivasan_MS_Thesis/Data/Data_Analysis_Results/OpenSwathWorkflow_3_results/Alphatims_Raw_IM_Extraction/filtered_diapasef_results_3_1.tsv")

