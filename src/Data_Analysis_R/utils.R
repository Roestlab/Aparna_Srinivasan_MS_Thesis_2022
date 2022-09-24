unimod_formatting <- function(peptide){
  peptide = gsub("_", "", peptide)
  peptide = gsub("\\(Acetyl \\(Protein N-term\\)\\)", "\\.\\(UniMod:1\\)", peptide)
  peptide = gsub("\\(Phospho \\(STY\\)\\)", "\\(UniMod:21\\)", peptide)
  peptide = gsub("\\(Oxidation \\(M\\)\\)", "\\(UniMod:35\\)", peptide)
  return(peptide)
}

codename_to_unimod_formatting <- function(peptide){
  peptide = gsub("\\(Acetylation\\)", "\\.\\(UniMod:1\\)", peptide)
  peptide = gsub("\\(Phospho\\)", "\\(UniMod:21\\)", peptide)
  peptide = gsub("\\(Oxidation\\)", "\\(UniMod:35\\)", peptide)
  peptide = gsub("\\(Carbamidomethyl\\)", "\\(UniMod:4\\)",peptide)
  return(peptide)
}

phospho_formatting <- function(peptide){
  peptide = gsub("\\.\\(UniMod:1\\)", "", peptide)
  peptide = gsub("\\(UniMod:35\\)", "", peptide)
  peptide = gsub("\\(UniMod:4\\)", "", peptide)
  return(peptide)
}

find_unique_ids <- function(df, pep_cutoff, peptide_name, pep_column){
  df %>% filter(get(pep_column) <= pep_cutoff) %>% group_by(method, mpool_gt) %>% 
    summarise(unique_ids = n_distinct(get(peptide_name))) %>% mutate(pep_cutoff = pep_cutoff) -> df_unique_obs
  
  return(df_unique_obs)
}

ipf_pep_pass <- function(ipf_pep_x, ipf_pep_y){
  if(ipf_pep_x <= 0.05 & ipf_pep_y <= 0.05){
    return(2)
  }
  else if(ipf_pep_x > 0.05 & ipf_pep_y > 0.05){
    return(0)
  }
  else if(ipf_pep_x <= 0.05 & ipf_pep_y > 0.05){
    return(1)
  }
  else if(ipf_pep_x > 0.05 & ipf_pep_y <= 0.05){
    return(1)
  }
}

ms2_pep_pass <- function(ms2_pep_x, ms2_pep_y){
  if(ms2_pep_x <= 0.05 & ms2_pep_y <= 0.05){
    return(2)
  }
  else if(ms2_pep_x > 0.05 & ms2_pep_y > 0.05){
    return(0)
  }
  else if(ms2_pep_x <= 0.05 & ms2_pep_y > 0.05){
    return(1)
  }
  else if(ms2_pep_x > 0.05 & ms2_pep_y <= 0.05){
    return(1)
  }
}


make_mcnemar_matrix <- function(df){
  mcnemar_matrix <- matrix(c(0,0,0,0), nrow=2, dimnames = list("Tims - On" = c("Identified", "Not Identified"), "Tims - Off" = c("Identified", "Not Identified")))
  mcnemar_matrix[1,1] = (df %>% filter(timson == 1 & timsoff == 1))$unique_phosphopeptide_ids
  mcnemar_matrix[1,2] = (df %>% filter(timson == 1 & timsoff == 0))$unique_phosphopeptide_ids
  mcnemar_matrix[2,1] = (df %>% filter(timson == 0 & timsoff == 1))$unique_phosphopeptide_ids
  mcnemar_matrix[2,2] = (df %>% filter(timson == 0 & timsoff == 0))$unique_phosphopeptide_ids
  return(mcnemar_matrix)
}

make_contingency_table <- function(lib_df, timsoff_df, timson_df, ramp){
  lib_df %>% select(Sequence, FullPeptideName, FullPeptideName_phospho) %>% distinct() -> p.data1
  
  timsoff_df %>% filter(mpool_gt == TRUE)  %>% filter(ms2_pep < 0.05) %>% filter(ipf_pep < 0.05) %>% select(FullPeptideName) %>% distinct() %>% mutate(timsoff = 1) -> p.data2
  
  timson_df %>% filter(ramp_time == ramp) %>% filter(mpool_gt == TRUE)  %>% filter(ms2_pep < 0.05) %>%  filter(ipf_pep < 0.05) %>% select(FullPeptideName) %>% distinct() %>% mutate(timson = 1) -> p.data3
  
  left_join(p.data1, p.data2, by="FullPeptideName") -> p.data1
  left_join(p.data1, p.data3, by="FullPeptideName") -> p.data1
  p.data1$timsoff %>% replace_na(0) -> p.data1$timsoff
  p.data1$timson %>% replace_na(0) -> p.data1$timson
  
  p.data1 <- p.data1 %>% 
    group_by(FullPeptideName_phospho) %>% summarise(sum(timsoff), sum(timson)) %>%
    mutate(timsoff = if_else(`sum(timsoff)` >= 1, 1,0 )) %>% 
    mutate(timson = if_else(`sum(timson)` >= 1, 1,0 )) %>% 
    select(-`sum(timsoff)`, -`sum(timson)`) %>% 
    group_by(timson, timsoff) %>% summarise(n_distinct(FullPeptideName_phospho)) %>% 
    rename(unique_phosphopeptide_ids = `n_distinct(FullPeptideName_phospho)`) %>% ungroup()
  
  return(p.data1)
}

make_phospho_contingency_table <- function(lib_df, timsoff_df, timson_df, ramp){
  lib_df %>% select(Sequence, FullPeptideName, FullPeptideName_phospho) %>% distinct() -> p.data1
  
  timsoff_df %>% filter(mpool_gt == TRUE)  %>% filter(ms2_pep < 0.05) %>% filter(ipf_pep < 0.05) %>% filter(grepl("UniMod:21", FullPeptideName)) %>% select(FullPeptideName) %>% distinct() %>% mutate(timsoff = 1) -> p.data2
  
  timson_df %>% filter(ramp_time == ramp) %>% filter(mpool_gt == TRUE)  %>% filter(ms2_pep < 0.05) %>%  filter(ipf_pep < 0.05) %>% filter(grepl("UniMod:21", FullPeptideName)) %>% select(FullPeptideName) %>% distinct() %>% mutate(timson = 1) -> p.data3
  
  left_join(p.data1, p.data2, by="FullPeptideName") -> p.data1
  left_join(p.data1, p.data3, by="FullPeptideName") -> p.data1
  p.data1$timsoff %>% replace_na(0) -> p.data1$timsoff
  p.data1$timson %>% replace_na(0) -> p.data1$timson
  
  p.data1 <- p.data1 %>% 
    group_by(FullPeptideName_phospho) %>% summarise(sum(timsoff), sum(timson)) %>%
    mutate(timsoff = if_else(`sum(timsoff)` >= 1, 1,0 )) %>% 
    mutate(timson = if_else(`sum(timson)` >= 1, 1,0 )) %>% 
    select(-`sum(timsoff)`, -`sum(timson)`) %>% 
    group_by(timson, timsoff) %>% summarise(n_distinct(FullPeptideName_phospho)) %>% 
    rename(unique_phosphopeptide_ids = `n_distinct(FullPeptideName_phospho)`) %>% ungroup()
  
  return(p.data1)
}

identification_table <- function(lib_df, timsoff_df, timson_df, ramp){
  lib_df %>% select(Sequence, FullPeptideName, FullPeptideName_phospho, pool) %>% distinct() -> p.data1
  
  timsoff_df %>% filter(mpool_gt == TRUE) %>% filter(ms2_pep < 0.05) %>% filter(ipf_pep < 0.05) %>% select(FullPeptideName) %>% distinct() %>% mutate(timsoff = 1) -> p.data2
  
  timson_df %>% filter(ramp_time == ramp) %>% filter(mpool_gt == TRUE) %>%  filter(ms2_pep < 0.05) %>% filter(ipf_pep < 0.05) %>% select(FullPeptideName) %>% distinct() %>% mutate(timson = 1) -> p.data3
  
  left_join(p.data1, p.data2, by="FullPeptideName") -> p.data1
  left_join(p.data1, p.data3, by="FullPeptideName") -> p.data1
  p.data1$timsoff %>% replace_na(0) -> p.data1$timsoff
  p.data1$timson %>% replace_na(0) -> p.data1$timson
  
  return(p.data1)
}

format_modified_sequence <- function(mod_seq){
  mod_seq = gsub("_", "", mod_seq)
  mod_seq = gsub("\\(Acetyl \\(Protein N-term\\)\\)", "", mod_seq)
  mod_seq = gsub("\\(Phospho \\(STY\\)\\)", "\\(UniMod:21\\)", mod_seq)
  mod_seq = gsub("\\(Oxidation \\(M\\)\\)", "", mod_seq)
  return(mod_seq)
}

phospho_position <- function(FullPeptideName_phospho){
  position = gregexpr("\\(", FullPeptideName_phospho)[[1]][1] -1 
  return(position)
}

histidine_position <- function(Sequence){
  position = gregexpr("H", Sequence)[[1]][1] 
  if(position == -1){
    return(0)
  }
  else{
    return(position)
  }
}

prolines_between_phosphosites <- function(Sequence, phosphosite_x, phosphosite_y){
  substring = substr(Sequence, pmin(phosphosite_x,phosphosite_y), pmax(phosphosite_x,phosphosite_y))
  n_p = str_count(substring, "P")
  return(n_p)
}

aa_next_to_phosphosite <- function(Sequence, phosphosite_x, phosphosite_y){
  substring = substr(Sequence, pmin(phosphosite_x,phosphosite_y), pmax(phosphosite_x,phosphosite_y))
  n_p = str_count(substring, "P")
  return(n_p)
}
