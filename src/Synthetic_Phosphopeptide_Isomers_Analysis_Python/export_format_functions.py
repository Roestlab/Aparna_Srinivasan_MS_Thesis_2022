import re
import numpy as np
import pandas as pd
from pathlib import Path, PurePosixPath 
import sys
import os
import sqlite3

def mapping_func(grouped_peptide_df):
    grouped_peptide_df = grouped_peptide_df[["ID", "MODIFIED_SEQUENCE"]]
    grouped_peptide_df["ID_TYPE"] = grouped_peptide_df["MODIFIED_SEQUENCE"].apply(lambda x: "UNIMOD_ID" if "UniMod" in x else "CODENAME_ID")
    grouped_peptide_df = grouped_peptide_df[["ID_TYPE", "ID"]]
    
    if len(grouped_peptide_df) == 1:
        if grouped_peptide_df["ID_TYPE"].item() == "CODENAME_ID":
            grouped_peptide_df = pd.concat([grouped_peptide_df, pd.DataFrame({"ID": [-1], "ID_TYPE":"UNIMOD_ID"})])
        else:
            grouped_peptide_df = pd.concat([grouped_peptide_df, pd.DataFrame({"ID": [-1], "ID_TYPE":"CODENAME_ID"})])
    
    grouped_peptide_df = pd.pivot_table(grouped_peptide_df, columns="ID_TYPE")
    return grouped_peptide_df

def make_mapping_table(osw_file):
    peptide_df = pd.read_sql_query("SELECT * FROM PEPTIDE", osw_file)
    peptide_df["tmp"] = peptide_df["MODIFIED_SEQUENCE"].replace("\\(\\w+\\)|\\(\\w+[:]\\d+\\)", "(@)", regex=True)
    peptide_df = peptide_df.groupby("tmp")
    peptide_df = peptide_df.apply(mapping_func) 
    peptide_df.columns.name = None
    peptide_df = peptide_df.reset_index(drop=True)

    return peptide_df

def add_mapping_table(osw):
    osw_file = sqlite3.connect(osw)
    df = make_mapping_table(osw_file)
    df.to_sql("UNIMOD_CODENAME_MAPPING", osw_file, index=False)
    osw_file.close()

def modification_only_phospho(modifiedsequence):
    """ only keep phosphorylation modification from modified sequence in unimod format"""
    modifiedsequence = modifiedsequence.replace("(UniMod:35)", "")
    modifiedsequence = modifiedsequence.replace("(UniMod:4)",  "")
    modifiedsequence = modifiedsequence.replace(".(UniMod:1)",  "")
    return modifiedsequence

def custom_pyprophet_export(osw):
    osw_file = sqlite3.connect(osw)

    export_query = '''
                    SELECT 
                        PRECURSOR.ID AS precursor_id,
                            PEPTIDE.ID AS peptide_id,
                            PEPTIDE.MODIFIED_SEQUENCE as FullPeptideName,
                            PEPTIDE.UNMODIFIED_SEQUENCE as Sequence,
                            PRECURSOR.PRECURSOR_MZ as precursor_mz,
                            RUN.FILENAME as filename,
                            PRECURSOR.CHARGE as precursor_charge,
                            FEATURE.ID as feature_id,
                            FEATURE.EXP_RT as RT,
                            FEATURE.LEFT_WIDTH as leftWidth,
                            FEATURE.RIGHT_WIDTH as rightWidth,
                            FEATURE.EXP_IM as IM,
                            FEATURE_MS2.AREA_INTENSITY as Intensity,
                            SCORE_MS2.RANK as peak_group_rank,
                            SCORE_MS2.QVALUE as ms2_m_score,
                            SCORE_MS2.PEP as ms2_pep,
                            SCORE_IPF.QVALUE as ipf_m_score,
                            SCORE_IPF.PEP as ipf_pep
                        FROM PRECURSOR
                        INNER JOIN PRECURSOR_PEPTIDE_MAPPING ON PRECURSOR.ID = PRECURSOR_PEPTIDE_MAPPING.PRECURSOR_ID
                        INNER JOIN PEPTIDE ON PRECURSOR_PEPTIDE_MAPPING.PEPTIDE_ID = PEPTIDE.ID
                        INNER JOIN FEATURE ON FEATURE.PRECURSOR_ID = PRECURSOR.ID
                        INNER JOIN RUN ON RUN.ID = FEATURE.RUN_ID
                        LEFT JOIN FEATURE_MS2 ON FEATURE_MS2.FEATURE_ID = FEATURE.ID
                        LEFT JOIN SCORE_MS2 ON SCORE_MS2.FEATURE_ID = FEATURE.ID
                        LEFT JOIN  (
                            SELECT 
                                SCORE_IPF.FEATURE_ID,
                                SCORE_IPF.PEPTIDE_ID,
                                UNIMOD_CODENAME_MAPPING.UNIMOD_ID,
                                SCORE_IPF.PRECURSOR_PEAKGROUP_PEP,
                                SCORE_IPF.QVALUE,
                                SCORE_IPF.PEP
                            FROM SCORE_IPF
                            INNER JOIN UNIMOD_CODENAME_MAPPING ON UNIMOD_CODENAME_MAPPING.CODENAME_ID = SCORE_IPF.PEPTIDE_ID
                        ) AS SCORE_IPF ON (SCORE_IPF.FEATURE_ID = FEATURE.ID AND SCORE_IPF.UNIMOD_ID = PEPTIDE.ID)
                        INNER JOIN PEPTIDE AS PEPTIDE_IPF ON SCORE_IPF.PEPTIDE_ID = PEPTIDE_IPF.ID
						INNER JOIN UNIMOD_CODENAME_MAPPING ON UNIMOD_CODENAME_MAPPING.UNIMOD_ID = PEPTIDE.ID
                        WHERE PRECURSOR.DECOY = 0
                        AND SCORE_MS2.RANK IS NOT NULL 
                    '''

    output = pd.read_sql_query(export_query, osw_file)
    return output

def get_pool_name(filename):
    return re.search("M(\\d{1,2})", filename).group(0)

def get_ramp_time(filename):
    if re.search("DIA-PaSEF", filename):
        if re.search("(?<=PaSEF)(\\d{3})(?=_)",filename) == None: 
            return 100.0
        else:
            return float(re.search("(?<=PaSEF)(\\d{3})(?=_)",filename).group(0)) 
    elif re.search("diaPASEF", filename):
        return float(re.search("(?<=\\_)(\\d{3})(?=ms)", filename).group(0))
    else:
        return float(re.search("(?<=\\-)(\\d{3})(?=ms)", filename).group(0))

def modification_codename_unimod(modifiedsequence):
    """Reformat the modified peptide sequence column in the MQ evidence.txt output for the probability extraction function"""
    modifiedsequence = modifiedsequence.replace("(Phospho)", "(UniMod:21)")
    modifiedsequence = modifiedsequence.replace("(Oxidation)", "(UniMod:35)")
    modifiedsequence = modifiedsequence.replace("(Carbamidomethyl)", "(UniMod:4)")
    modifiedsequence = modifiedsequence.replace("(Acetyl)", "(UniMod:1)")

    return modifiedsequence

def remove_modifications(FullPeptideName):
    sequence = FullPeptideName.replace("(UniMod:35)", "")
    sequence = sequence.replace("(UniMod:4)",  "")
    sequence = sequence.replace(".(UniMod:1)",  "")
    sequence = sequence.replace("(UniMod:21)",  "")
    return sequence

def fragments_with_phosphorylation(FullPeptideName):
    sequence = remove_modifications(FullPeptideName)
    l = len(sequence)
    FullPeptideName_phospho = modification_only_phospho(FullPeptideName)

    n = FullPeptideName_phospho.find("(")

    fragment_list = []
    for i in range(n,l):
        fragment_list.append(f"b{i}")

    for i in range(l-n+1, l):
        fragment_list.append(f"y{i}")
    
    return fragment_list

def phosphorylation_position(FullPeptideName):
    FullPeptideName_phospho = modification_only_phospho(FullPeptideName)
    n = FullPeptideName_phospho.find("(")
    return n
