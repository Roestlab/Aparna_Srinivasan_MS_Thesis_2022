library(tidyverse)
library(data.t)
source("/Users/aparnasrinivasan/Desktop/MScThesis/Github/Synthetic_phosphopeptide_isomers_TimsTOF/Analysis/utils.R")

ground_truth <- data.table::fread("/Users/aparnasrinivasan/Desktop/MScThesis/Srinivasan_MS_Thesis/Data/ground_truth_table.tsv")
ground_truth$p_res <- mapply(phosphorylatable_residues, ground_truth$Sequence)

evidence_annotated <- data.table::fread("/Users/aparnasrinivasan/Desktop/MScThesis/Srinivasan_MS_Thesis/Data/Data_Analysis_Results/Upool_spectral_library_generation/evidence_annotated.tsv") %>% 
  left_join(ground_truth %>% select(Sequence) %>% mutate(protein_in_gt = TRUE) %>% rename(Proteins = Sequence), by="Proteins") %>% 
  filter(!grepl("REMOVE", Proteins)) %>% 
  mutate(FullPeptideName = unimod_formatting(`Modified sequence`))

isomer_pair_df <- data.table::fread("/Users/aparnasrinivasan/Desktop/MScThesis/Srinivasan_MS_Thesis/Data/Data_Analysis_Results/OpenSwathWorkflow_3_results/Isomers_Analysis/co_eluting_isomer_pairs.tsv") #%>% 
 #  mutate(delta_IM_transition = abs(transition_IM.x- transition_IM.y))# %>% 
 # filter(intefering_features==2)

#optional
isomer_pairs <- c("1388467681012992100 2913950545086492932","541526734563420968 8620910759637318001", "754444565133828392 8002785594889782764")
isomer_pair_df %>% filter(isomer_pair %in% isomer_pairs) -> isomer_pair_df


df <- evidence_annotated %>% filter(ground_truth_upool==TRUE) %>% filter(SiteProb==1) %>% 
  select(Sequence, FullPeptideName, `Raw file`, Experiment, Charge, `m/z`, `Retention time`, `1/K0`, `PEP`) %>% 
  unique() %>% group_by(`Raw file`, FullPeptideName, Charge) %>% 
  slice(which.min(PEP)) %>% ungroup() %>% 
  rename(filename = `Raw file`, mz = `m/z`, RT = `Retention time`, IM = `1/K0`)  %>% 
  mutate(RT = RT*60)

df.1 <- isomer_pair_df %>% select(feature_id.x, feature_id.y, FullPeptideName.x, FullPeptideName.y, library_IM.x, library_IM.y, Charge, mz) %>% 
  mutate(FullPeptideName.x = gsub("\\(UniMod:4\\)", "", FullPeptideName.x), FullPeptideName.y = gsub("\\(UniMod:4\\)", "", FullPeptideName.y))  # max quant doesn't explicitly report carb modifications on cysteine 




df.1 %>% left_join(df %>% select(FullPeptideName, Charge, filename, RT, IM) %>% 
          rename(FullPeptideName.x = FullPeptideName, filename.x = filename, RT.x = RT, IM.x = IM), by=c("FullPeptideName.x", "Charge")) %>%  
          left_join(df %>% select(FullPeptideName, Charge,  filename, RT, IM) %>% 
                      rename(FullPeptideName.y = FullPeptideName, filename.y = filename, RT.y = RT, IM.y = IM), by=c("FullPeptideName.y", "Charge")) %>% 
   mutate(Upool_pair = paste(pmin(filename.x, filename.y), pmax(filename.x, filename.y))) %>% 
  mutate(FullPeptideName.x = gsub("C", "C\\(UniMod:4\\)", FullPeptideName.x), FullPeptideName.y = gsub("C", "C\\(UniMod:4\\)", FullPeptideName.y)) -> df.1  # add back carb modification 




write_tsv(df.1, "/Users/aparnasrinivasan/Desktop/MScThesis/Srinivasan_MS_Thesis/Data/Data_Analysis_Results/OpenSwathWorkflow_3_results/Thesis_figures/DDA_mobilogram_overlay_dataframe.tsv")

