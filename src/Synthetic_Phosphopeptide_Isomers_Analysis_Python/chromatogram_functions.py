import sqlite3
import PyMSNumpress as pyms
import pandas as pd
import zlib
from scipy.signal import savgol_filter
import matplotlib.pyplot as plt
import seaborn as sns
from export_format_functions import *

# According to OpenMS documentation ,
# sqmass data types: 0 =mz, 1 = intensity, 2 = retention time 
# sqmass compression types 0 = no compression, 1= zlib, 2 = np-linear, 3 = np-slof, 4 = np-pic, 5 = np-linear + zlib
# 6 = np-slof + zlib, 7 = np-pic + zlib (as described in the MSNumpress paper) 
# intensity compression is 6, rt compression is 5 
# pandas > 1.2 is not supported in python 3.6
# however i must use python3.6 in order to use MSNumpress 
# score pyprophet files already have the unimod codename mapping 


#plt.figure(figsize=(12,8), dpi=300) 
# might wanna add this to fix the plot to legend ratio 

def export_chromatogram_df(sqmass_file):
    sqmass_file = sqlite3.connect(sqmass_file)

    sqmass_export_query  = '''
                        SELECT 
                        DATA.CHROMATOGRAM_ID AS chromatogram_id,
                        DATA.COMPRESSION AS compression_type,
                        DATA.DATA_TYPE AS data,
                        DATA.DATA AS array, 
                            CHROMATOGRAM.NATIVE_ID AS transition_id
                        FROM DATA 
                        INNER JOIN CHROMATOGRAM ON CHROMATOGRAM.ID = DATA.CHROMATOGRAM_ID
                        '''
    chromatogram_data = pd.read_sql_query(sqmass_export_query, sqmass_file)# dtype={"transition_id":int, "compression_type": int, "chromatogram_id": int, "data":int}) 
    chromatogram_data = chromatogram_data.astype({"transition_id":int}, copy=False)
    sqmass_file.close()
    return chromatogram_data

# the following ony export identification transitions 
def export_identification_transition_df(osw_file):
    osw_file = sqlite3.connect(osw_file)
    
    osw_export_query = '''
                   SELECT
                   FEATURE.ID as feature_id,
                        FEATURE.EXP_RT as feature_retention_time,
                        FEATURE.EXP_IM as feature_IM, 
                        PEPTIDE.MODIFIED_SEQUENCE as FullPeptideName,
                        PRECURSOR.PRECURSOR_MZ as precursor_mz, 
                        TRANSITION_PRECURSOR_MAPPING.TRANSITION_ID as transition_id_precursor,
                        TRANSITION.PRODUCT_MZ as transition_mz,
                        TRANSITION.ANNOTATION as transition_annotation,
                        TRANSITION.DETECTING as transition_detecting,
                        TRANSITION.IDENTIFYING as transition_identifying
                   FROM FEATURE 
                   INNER JOIN PRECURSOR_PEPTIDE_MAPPING ON PRECURSOR_PEPTIDE_MAPPING.PRECURSOR_ID = FEATURE.PRECURSOR_ID
                   INNER JOIN UNIMOD_CODENAME_MAPPING ON UNIMOD_CODENAME_MAPPING.UNIMOD_ID = PRECURSOR_PEPTIDE_MAPPING.PEPTIDE_ID
                   INNER JOIN TRANSITION_PRECURSOR_MAPPING ON TRANSITION_PRECURSOR_MAPPING.PRECURSOR_ID = FEATURE.PRECURSOR_ID
                   INNER JOIN TRANSITION_PEPTIDE_MAPPING ON TRANSITION_PEPTIDE_MAPPING.PEPTIDE_ID = UNIMOD_CODENAME_MAPPING.CODENAME_ID AND TRANSITION_PEPTIDE_MAPPING.TRANSITION_ID = transition_id_precursor
                   INNER JOIN TRANSITION ON TRANSITION_PEPTIDE_MAPPING.TRANSITION_ID = TRANSITION.ID
                   INNER JOIN PRECURSOR ON FEATURE.PRECURSOR_ID = PRECURSOR.ID
                   INNER JOIN PEPTIDE ON PEPTIDE.ID = UNIMOD_CODENAME_MAPPING.CODENAME_ID
                   WHERE TRANSITION.DECOY = 0 AND PRECURSOR.DECOY = 0 
                   '''
    
    feature_df = pd.read_sql_query(osw_export_query, osw_file) # dtype={"transition_id":int, "feature_id":int, "precursor_id":int, "transition_mz":float, "transition_annotation":str, "transition_detecting":int, "transition_identifying":int})
    feature_df["FullPeptideName"] = feature_df["FullPeptideName"].apply(lambda x: modification_codename_unimod(x))
    feature_df = feature_df.rename(columns={"transition_id_precursor": "transition_id"})
    osw_file.close()
    return feature_df

def export_detection_transition_df(osw_file):
    osw_file = sqlite3.connect(osw_file)
    
    osw_export_query = '''
                   SELECT
                   FEATURE.ID as feature_id,
                        FEATURE.EXP_RT as feature_retention_time,
                        FEATURE.EXP_IM as feature_IM, 
                        PEPTIDE.MODIFIED_SEQUENCE as FullPeptideName, 
                        PRECURSOR.PRECURSOR_MZ as precursor_mz, 
                        TRANSITION_PRECURSOR_MAPPING.TRANSITION_ID as transition_id_precursor,
                        TRANSITION.PRODUCT_MZ as transition_mz,
                        TRANSITION.ANNOTATION as transition_annotation,
                        TRANSITION.DETECTING as transition_detecting,
                        TRANSITION.IDENTIFYING as transition_identifying
                   FROM FEATURE 
                   INNER JOIN PRECURSOR_PEPTIDE_MAPPING ON PRECURSOR_PEPTIDE_MAPPING.PRECURSOR_ID = FEATURE.PRECURSOR_ID
                   INNER JOIN UNIMOD_CODENAME_MAPPING ON UNIMOD_CODENAME_MAPPING.UNIMOD_ID = PRECURSOR_PEPTIDE_MAPPING.PEPTIDE_ID
                   INNER JOIN TRANSITION_PRECURSOR_MAPPING ON TRANSITION_PRECURSOR_MAPPING.PRECURSOR_ID = FEATURE.PRECURSOR_ID
                   INNER JOIN TRANSITION ON TRANSITION.ID = TRANSITION_PRECURSOR_MAPPING.TRANSITION_ID
                   INNER JOIN PRECURSOR ON FEATURE.PRECURSOR_ID = PRECURSOR.ID
                   INNER JOIN PEPTIDE ON PEPTIDE.ID = UNIMOD_CODENAME_MAPPING.CODENAME_ID
                   WHERE TRANSITION.DECOY = 0 AND PRECURSOR.DECOY = 0 AND TRANSITION.DETECTING = 1
                   '''
    
    feature_df = pd.read_sql_query(osw_export_query, osw_file) # dtype={"transition_id":int, "feature_id":int, "precursor_id":int, "transition_mz":float, "transition_annotation":str, "transition_detecting":int, "transition_identifying":int})
    feature_df["FullPeptideName"] = feature_df["FullPeptideName"].apply(lambda x: modification_codename_unimod(x))
    feature_df = feature_df.rename(columns={"transition_id_precursor": "transition_id"})
    osw_file.close()
    return feature_df
  
def export_all_transitions(osw_file):
    
    df_1 = export_detection_transition_df(osw_file)
    df_2 = export_identification_transition_df(osw_file)

    transition_df = pd.concat([df_1,df_2], axis=0)

    return transition_df 

def peptide_df(feature_id, feature_df, chromatogram_df):
    feature = feature_df.loc[feature_df["feature_id"]==feature_id]
    feature = feature.loc[feature.groupby("transition_annotation")["transition_identifying"].idxmax()]

    chromatogram_data = chromatogram_df.loc[chromatogram_df["transition_id"].isin(feature["transition_id"])]
    peptide_data = feature.merge(chromatogram_data, on="transition_id")
    return peptide_data


def decode_intensity(intensity):
    intensity = zlib.decompress(intensity)
    decode_intensity = []
    pyms.decodeSlof(intensity, decode_intensity)
    return decode_intensity

def decode_retention_time(retention_time):
    retention_time = zlib.decompress(retention_time)
    decode_rt = []
    pyms.decodeLinear(retention_time, decode_rt)
    return decode_rt


def plot_peptide_chromatogram(feature, feature_df, chrom_df):
    
    chromatogram_df = peptide_df(feature, feature_df, chrom_df)
    peptide = chromatogram_df["FullPeptideName"].unique().item()

    feature_rt = chromatogram_df["feature_retention_time"].unique().item()

    clrs_1 = sns.color_palette('Blues', n_colors=len(set(chromatogram_df.loc[chromatogram_df["transition_detecting"] ==1]["transition_id"])))
    clrs_2 = sns.color_palette('Greens', n_colors=len(set(chromatogram_df.loc[chromatogram_df["transition_identifying"] ==1]["transition_id"])))

    i=0
    j=0
    
    plt.figure(figsize=(10,8))
    plt.style.use("seaborn")

    for transition in set(chromatogram_df["transition_id"]):
            intensity = chromatogram_df[(chromatogram_df["transition_id"]==transition) & (chromatogram_df["data"]==1)]["array"].item()
            retention_time = chromatogram_df[(chromatogram_df["transition_id"]==transition) & (chromatogram_df["data"]==2)]["array"].item()
            transition_annotation = chromatogram_df[(chromatogram_df["transition_id"]==transition) & (chromatogram_df["data"]==2)]["transition_annotation"].item()
            decoded_rt = decode_retention_time(retention_time)
            decoded_intensity = decode_intensity(intensity)
            smoothed_intensity = savgol_filter(decoded_intensity, 11,3)

            
            if chromatogram_df[(chromatogram_df["transition_id"]==transition) & (chromatogram_df["data"]==1)]["transition_detecting"].item() ==1:
                plt.plot(decoded_rt, smoothed_intensity, color=clrs_1[i], linestyle="--", alpha=1, label=f"{transition_annotation}")
                i +=1 
            else: 
                plt.plot(decoded_rt, smoothed_intensity, color=clrs_2[j], linestyle="--", alpha=1, label=f"{transition_annotation}")
                j +=1 

    
    plt.title(f"{peptide}", loc="left")
    plt.xlabel("Retention time (s)")
    plt.ylabel("Intensity")
    plt.legend(bbox_to_anchor=(1.05,1), loc="upper left")
    plt.axhline(y=0, color = "black", linewidth=1)
    plt.axvline(x=feature_rt, color="red", linewidth=0.5)
    plt.tight_layout()
    plt.savefig(f"{feature}_chromatogram.png", dpi=300)
    plt.clf()

def plot_unique_transitions(feature_1, feature_2, feature_df, chrom_df):
    
    chromatogram_df = pd.concat([peptide_df(feature_1, feature_df, chrom_df), peptide_df(feature_2, feature_df, chrom_df)])
    chromatogram_df["n_features"] = (chromatogram_df.groupby(["transition_mz", "transition_annotation"])["feature_id"].transform("nunique"))
    chromatogram_df = chromatogram_df.loc[chromatogram_df["n_features"] != 2]
    
    clrs_1 = sns.color_palette('Blues', n_colors = len(set(chromatogram_df[chromatogram_df["feature_id"]==feature_1]["transition_id"])))
    clrs_2 = sns.color_palette('Greens', n_colors = len(set(chromatogram_df[chromatogram_df["feature_id"]==feature_2]["transition_id"])))
    
    peptides = chromatogram_df["FullPeptideName"].unique()
    i = 0 
    j = 0

    plt.figure(figsize=(10,8))
    plt.style.use("seaborn")

    for transition in set(chromatogram_df["transition_id"]):
            intensity = chromatogram_df[(chromatogram_df["transition_id"]==transition) & (chromatogram_df["data"]==1)]["array"].item()
            retention_time = chromatogram_df[(chromatogram_df["transition_id"]==transition) & (chromatogram_df["data"]==2)]["array"].item()
            peptide = chromatogram_df[(chromatogram_df["transition_id"]==transition) & (chromatogram_df["data"]==2)]["FullPeptideName"].item()
            transition_annotation = chromatogram_df[(chromatogram_df["transition_id"]==transition) & (chromatogram_df["data"]==2)]["transition_annotation"].item()
            decoded_rt = decode_retention_time(retention_time)
            decoded_intensity = decode_intensity(intensity)
            smoothed_intensity = savgol_filter(decoded_intensity, 11,3)
           
            if peptide == peptides[0]: 
                plt.plot(decoded_rt, smoothed_intensity, color=clrs_1[i], linestyle="--", alpha=1, label=f"{transition_annotation}_{peptide}")
                i+=1 
            else: 
                plt.plot(decoded_rt, smoothed_intensity, color=clrs_2[j], linestyle="--", alpha=1, label=f"{transition_annotation}_{peptide}")
                j+=1 
    
    plt.title(f"{peptides[0]} and {peptides[1]}", loc="left")
    plt.xlabel("Retention time (s)")
    plt.ylabel("Intensity")
    plt.legend(bbox_to_anchor=(1.05,1), loc="upper left")
    plt.axhline(y=0, color = "black", linewidth=1)
    plt.tight_layout()
    plt.savefig(f"{feature_1}_{feature_2}_unique_fragments_chromatogram.png", dpi=300)
    plt.clf()
