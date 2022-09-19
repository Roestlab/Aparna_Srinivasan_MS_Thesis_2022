import sqlite3
import PyMSNumpress as pyms
import pandas as pd
import zlib
from scipy.signal import savgol_filter
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from export_format_functions import *
import seaborn as sns

# According to OpenMS documentation ,
# sqmass data types: 0 =mz, 1 = intensity, 2 = retention time 
# sqmass compression types 0 = no compression, 1= zlib, 2 = np-linear, 3 = np-slof, 4 = np-pic, 5 = np-linear + zlib
# 6 = np-slof + zlib, 7 = np-pic + zlib (as described in the MSNumpress paper) 
# intensity compression is 6, rt compression is 5 
# pandas > 1.2 is not supported in python 3.6
# however i must use python3.6 in order to use MSNumpress so...
# score pyprophet files already have the unimod codename mapping 

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
    #chromatogram_data = chromatogram_data.astype({"transition_id":int}, copy=False)
    sqmass_file.close()
    return chromatogram_data

def export_transition_df(osw_file):
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

def export_transition_chromatogram_df(sqmass_file):
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
    chromatogram_data = chromatogram_data.loc[chromatogram_data["transition_id"].str.contains("Precursor")==False]
    chromatogram_data = chromatogram_data.astype({"transition_id":int}, copy=False)
    sqmass_file.close()
    return chromatogram_data


def export_precursor_chromatogram_df(sqmass_file):
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
    chromatogram_data = pd.read_sql_query(sqmass_export_query, sqmass_file) # dtype={"transition_id":int, "compression_type": int, "chromatogram_id": int, "data":int}) 
    sqmass_file.close()
    # keep rows that contain "precursor" in the transition id column 
    chromatogram_data = chromatogram_data.loc[chromatogram_data["transition_id"].str.contains("Precursor")==True]
    chromatogram_data["precursor_id"] = chromatogram_data["transition_id"].apply(lambda x: x.split("_")[0])
    chromatogram_data["Isotope"] = chromatogram_data["transition_id"].apply(lambda x: x.split("_")[2])
    chromatogram_data = chromatogram_data.astype({"precursor_id":int}, copy=False)
    return chromatogram_data

def peptide_df(feature_id, feature_df, chromatogram_df):
    feature = feature_df.loc[feature_df["feature_id"]==feature_id]
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


def merged_transition_df(feature_df, chrom_transition):
    chromatogram_transition_df = feature_df.merge(chrom_transition, on="transition_id") 
    return chromatogram_transition_df

def merged_precursor_df(feature_df, chrom_precursor):
    chromatogram_precursor_df = feature_df.drop(columns={"transition_id", "transition_mz", "transition_annotation", "transition_detecting", "transition_identifying", "transition_charge", "transition_annotation_nocharge"}, axis=1)
    chromatogram_precursor_df = chromatogram_precursor_df.drop_duplicates()
    chromatogram_precursor_df = chromatogram_precursor_df.merge(chrom_precursor, on="precursor_id")
    chromatogram_precursor_df = chromatogram_precursor_df.loc[chromatogram_precursor_df["Isotope"]=="i0"]
    return chromatogram_precursor_df



def plot_peptide_chromatogram(feature, feature_df, chrom_df):
    chromatogram_df = peptide_df(feature, feature_df, chrom_df)
    peptide = chromatogram_df["FullPeptideName"].unique().item()

    fig_handles = []
    colors = ["dodgerblue", "fuchsia"]
    transition_type = ["Identifying transition", "Detecting transition"]

    for i in range(2):
        fig_handles.append(mpatches.Patch(color=colors[i], label=f"{transition_type[i]}"))


    for transition in set(chromatogram_df["transition_id"]):
            intensity = chromatogram_df[(chromatogram_df["transition_id"]==transition) & (chromatogram_df["data"]==1)]["array"].item()
            retention_time = chromatogram_df[(chromatogram_df["transition_id"]==transition) & (chromatogram_df["data"]==2)]["array"].item()
            
            decoded_rt = decode_retention_time(retention_time)
            decoded_intensity = decode_intensity(intensity)
            smoothed_intensity = savgol_filter(decoded_intensity, 11,3)
           
            if chromatogram_df[(chromatogram_df["transition_id"]==transition) & (chromatogram_df["data"]==2)]["transition_identifying"].item() ==1:
                plt.plot(decoded_rt, smoothed_intensity, color=colors[0], linestyle="--", alpha=1)
            else: 
                plt.plot(decoded_rt, smoothed_intensity, color=colors[1], linestyle="--", alpha=1)
    

    plt.legend(handles=fig_handles)
    plt.title(f"{peptide}")
    plt.savefig(f"{feature}_chromatogram.png")
    plt.clf()
    

def plot_unique_transitions(feature_1, feature_2, feature_df, chrom_df):
    
    chromatogram_df = pd.concat([peptide_df(feature_1, feature_df, chrom_df), peptide_df(feature_2, feature_df, chrom_df)])
    chromatogram_df["n_features"] = (chromatogram_df.groupby(["transition_mz", "transition_annotation"])["feature_id"].transform("nunique"))
    chromatogram_df = chromatogram_df.loc[chromatogram_df["n_features"] != 2]
    
    fig_handles = []
    colors = ["dodgerblue", "fuchsia"]
    peptides = chromatogram_df["FullPeptideName"].unique()

    for i in range(2):
        fig_handles.append(mpatches.Patch(color=colors[i], label=f"{peptides[i]}"))


    for transition in set(chromatogram_df["transition_id"]):
            intensity = chromatogram_df[(chromatogram_df["transition_id"]==transition) & (chromatogram_df["data"]==1)]["array"].item()
            retention_time = chromatogram_df[(chromatogram_df["transition_id"]==transition) & (chromatogram_df["data"]==2)]["array"].item()
            peptide = chromatogram_df[(chromatogram_df["transition_id"]==transition) & (chromatogram_df["data"]==2)]["FullPeptideName"].item()

            decoded_rt = decode_retention_time(retention_time)
            decoded_intensity = decode_intensity(intensity)
            smoothed_intensity = savgol_filter(decoded_intensity, 11,3)
           
            if peptide == peptides[0]: 
                plt.plot(decoded_rt, smoothed_intensity, color=colors[0], linestyle="--", alpha=1)
            else: 
                plt.plot(decoded_rt, smoothed_intensity, color=colors[1], linestyle="--", alpha=1)
    

    plt.legend(handles=fig_handles)
    plt.savefig(f"{feature_1}_{feature_2}_unique_fragments_chromatogram.png")
    plt.clf()


def plot_all_transitions(feature_1, feature_2, feature_df, chrom_df):
    chromatogram_df = pd.concat([peptide_df(feature_1, feature_df, chrom_df), peptide_df(feature_2, feature_df, chrom_df)])
    fig_handles = []
    colors = ["dodgerblue", "fuchsia"]
    peptides = chromatogram_df["FullPeptideName"].unique()

    for i in range(2):
        fig_handles.append(mpatches.Patch(color=colors[i], label=f"{peptides[i]}"))

    for transition in set(chromatogram_df["transition_id"]):
            intensity = chromatogram_df[(chromatogram_df["transition_id"]==transition) & (chromatogram_df["data"]==1)]["array"].item()
            retention_time = chromatogram_df[(chromatogram_df["transition_id"]==transition) & (chromatogram_df["data"]==2)]["array"].item()
            peptide = chromatogram_df[(chromatogram_df["transition_id"]==transition) & (chromatogram_df["data"]==2)]["FullPeptideName"].item()

            decoded_rt = decode_retention_time(retention_time)
            decoded_intensity = decode_intensity(intensity)
            smoothed_intensity = savgol_filter(decoded_intensity, 11,3)
           
            if peptide == peptides[0]: 
                plt.plot(decoded_rt, smoothed_intensity, color=colors[0], linestyle="--", alpha=1)
            else: 
                plt.plot(decoded_rt, smoothed_intensity, color=colors[1], linestyle="--", alpha=1)
    

    plt.legend(handles=fig_handles)
    plt.savefig(f"{feature_1}_{feature_2}_unique_fragments_chromatogram.png")
    plt.clf()

def export_identification_transition_df(osw_file):
    osw_file = sqlite3.connect(osw_file)
    osw_export_query = '''
                   SELECT
                   FEATURE.ID as feature_id,
                        FEATURE.EXP_RT as feature_retention_time,
                        FEATURE.EXP_IM as feature_IM, 
                        FEATURE.PRECURSOR_ID as precursor_id,
                        PEPTIDE.MODIFIED_SEQUENCE as FullPeptideName,
                        PRECURSOR.PRECURSOR_MZ as precursor_mz, 
                        PRECURSOR.LIBRARY_DRIFT_TIME as library_IM, 
                        TRANSITION_PRECURSOR_MAPPING.TRANSITION_ID as transition_id_precursor,
                        TRANSITION.PRODUCT_MZ as transition_mz,
                        TRANSITION.ANNOTATION as transition_annotation,
                        TRANSITION.DETECTING as transition_detecting,
                        TRANSITION.IDENTIFYING as transition_identifying,
                        TRANSITION.CHARGE as transition_charge
                   FROM FEATURE 
                   INNER JOIN PRECURSOR_PEPTIDE_MAPPING ON PRECURSOR_PEPTIDE_MAPPING.PRECURSOR_ID = FEATURE.PRECURSOR_ID
                   INNER JOIN UNIMOD_CODENAME_MAPPING ON UNIMOD_CODENAME_MAPPING.UNIMOD_ID = PRECURSOR_PEPTIDE_MAPPING.PEPTIDE_ID
                   INNER JOIN TRANSITION_PRECURSOR_MAPPING ON TRANSITION_PRECURSOR_MAPPING.PRECURSOR_ID = FEATURE.PRECURSOR_ID
                   INNER JOIN TRANSITION_PEPTIDE_MAPPING ON TRANSITION_PEPTIDE_MAPPING.PEPTIDE_ID = UNIMOD_CODENAME_MAPPING.CODENAME_ID AND TRANSITION_PEPTIDE_MAPPING.TRANSITION_ID = transition_id_precursor
                   INNER JOIN TRANSITION ON TRANSITION_PEPTIDE_MAPPING.TRANSITION_ID = TRANSITION.ID
                   INNER JOIN PRECURSOR ON FEATURE.PRECURSOR_ID = PRECURSOR.ID
                   INNER JOIN PEPTIDE ON PEPTIDE.ID = UNIMOD_CODENAME_MAPPING.CODENAME_ID
                   WHERE TRANSITION.DECOY = 0 AND PRECURSOR.DECOY = 0 AND TRANSITION.IDENTIFYING = 1
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
                        FEATURE.PRECURSOR_ID as precursor_id,
                        PEPTIDE.MODIFIED_SEQUENCE as FullPeptideName, 
                        PRECURSOR.PRECURSOR_MZ as precursor_mz, 
                        PRECURSOR.LIBRARY_DRIFT_TIME as library_IM, 
                        TRANSITION_PRECURSOR_MAPPING.TRANSITION_ID as transition_id_precursor,
                        TRANSITION.PRODUCT_MZ as transition_mz,
                        TRANSITION.ANNOTATION as transition_annotation,
                        TRANSITION.DETECTING as transition_detecting,
                        TRANSITION.IDENTIFYING as transition_identifying,
                        TRANSITION.CHARGE as transition_charge
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
    transition_df["transition_annotation_nocharge"] = transition_df["transition_annotation"].apply(lambda x: x[:-2])
    return transition_df 
 
def chromatogram_overlay_precursors_unique_phosphorylated_transitions_annotated_nosmoothing(feature_1, feature_2, chromatogram_transitions, chromatogram_precursors):

    peptide_df_1_transitions = chromatogram_transitions.loc[chromatogram_transitions["feature_id"]==feature_1]
    peptide_1 = peptide_df_1_transitions[peptide_df_1_transitions["feature_id"]==feature_1]["FullPeptideName"].unique().item()
    peptide_df_1_transitions = peptide_df_1_transitions.loc[peptide_df_1_transitions["transition_annotation_nocharge"].isin(fragments_with_phosphorylation(peptide_1))]

    peptide_df_2_transitions = chromatogram_transitions.loc[chromatogram_transitions["feature_id"]==feature_2]
    peptide_2 = peptide_df_2_transitions[peptide_df_2_transitions["feature_id"]==feature_2]["FullPeptideName"].unique().item()
    peptide_df_2_transitions = peptide_df_2_transitions.loc[peptide_df_2_transitions["transition_annotation_nocharge"].isin(fragments_with_phosphorylation(peptide_2))]

    peptide_df_transitions = pd.concat([peptide_df_1_transitions, peptide_df_2_transitions], axis=0)
    peptide_df_transitions = peptide_df_transitions.loc[peptide_df_transitions["transition_identifying"] ==1]
    peptide_df_transitions["n_features"] = (peptide_df_transitions.groupby(["transition_annotation_nocharge"])["feature_id"].transform("nunique")) 
    peptide_df_transitions = peptide_df_transitions.loc[peptide_df_transitions["n_features"] != 2] 
    
    peptide_df_precursors = chromatogram_precursors.loc[chromatogram_precursors["feature_id"].isin([feature_1, feature_2])]
    RT_1 = peptide_df_precursors[peptide_df_precursors["feature_id"]==feature_1]["feature_retention_time"].unique().item()
    RT_2 = peptide_df_precursors[peptide_df_precursors["feature_id"]==feature_2]["feature_retention_time"].unique().item()

    plt.figure(figsize=(10,8))
    plt.style.use("seaborn")

    plot_lines_1 = []
    plot_lines_2 = []
    feature_label_1 = []
    feature_label_2 = []

    feature_list = np.array([feature_1, feature_2])
    phosphorylation_positions = np.array([phosphorylation_position(peptide_1), phosphorylation_position(peptide_2)])
    phosphorylation_positions = np.argsort(phosphorylation_positions)
    features_sorted = feature_list[phosphorylation_positions]
    features_sorted = np.flip(features_sorted)


    for j in range(2):
        feature = features_sorted[j]
        peptide = peptide_df_transitions.loc[peptide_df_transitions["feature_id"]==feature]["FullPeptideName"].unique().item()
        transition_df_peptide = peptide_df_transitions.loc[peptide_df_transitions["feature_id"] == feature] 
        RT = peptide_df_transitions.loc[peptide_df_transitions["feature_id"] == feature]["feature_retention_time"].unique().item()
        if j ==1: 
            clrs = sns.color_palette('Blues', n_colors= len(set(transition_df_peptide["transition_id"]))+1)
            rt_color = "deepskyblue"
            plt.axvline(x=RT, color=rt_color, linewidth=1.5, linestyle=":")

            precursor_intensities = peptide_df_precursors[(peptide_df_precursors["feature_id"]==feature) & (peptide_df_precursors["data"]==1)]["array"].item()
            precursor_rt = peptide_df_precursors[(peptide_df_precursors["feature_id"]==feature) & (peptide_df_precursors["data"]==2)]["array"].item()
            decoded_precursor_intensities = decode_intensity(precursor_intensities)
            decoded_precursor_rt = decode_retention_time(precursor_rt)
            l0, = plt.plot(decoded_precursor_rt, decoded_precursor_intensities, color="black")
            plot_lines_2.append(l0)
            feature_label_2.append("precursor")

            i=0
            for transition in set(transition_df_peptide["transition_id"]):
                if len(transition_df_peptide) > 0:
                    transition_annotation = transition_df_peptide[transition_df_peptide["transition_id"]==transition]["transition_annotation"].unique().item()

                    transition_intensities = transition_df_peptide[(transition_df_peptide["transition_id"]==transition) & (transition_df_peptide["data"]==1)]["array"].item()
                    transition_rt = transition_df_peptide[(transition_df_peptide["transition_id"]==transition) & (transition_df_peptide["data"]==2)]["array"].item()

                    decoded_transition_intensities = decode_intensity(transition_intensities)
                    decoded_transition_rt = decode_retention_time(transition_rt)

                    l2, = plt.plot(decoded_transition_rt, decoded_transition_intensities, color=clrs[i+1], linestyle="-", alpha=1)
                    plot_lines_2.append(l2)
                    feature_label_2.append(str(transition_annotation))
                    i+=1
                else:
                    break 
            
            legend_2 = plt.legend(plot_lines_2, feature_label_2, title=f"{peptide}", bbox_to_anchor=(1.05,0.20), loc="center left",fontsize=7)

        
        else:
            clrs = sns.color_palette('Greens', n_colors= len(set(transition_df_peptide["transition_id"]))+1)
            rt_color = "forestgreen"
            plt.axvline(x=RT, color=rt_color, linewidth=1.5, linestyle=":")
            i=0
            
            precursor_intensities = peptide_df_precursors[(peptide_df_precursors["feature_id"]==feature) & (peptide_df_precursors["data"]==1)]["array"].item()
            precursor_rt = peptide_df_precursors[(peptide_df_precursors["feature_id"]==feature) & (peptide_df_precursors["data"]==2)]["array"].item()
            decoded_precursor_intensities = decode_intensity(precursor_intensities)
            decoded_precursor_rt = decode_retention_time(precursor_rt)
            l0, = plt.plot(decoded_precursor_rt, decoded_precursor_intensities, color="black")
            plot_lines_1.append(l0)
            feature_label_1.append("precursor")



            for transition in set(transition_df_peptide["transition_id"]):
                if len(transition_df_peptide) > 0:
                    transition_annotation = transition_df_peptide[transition_df_peptide["transition_id"]==transition]["transition_annotation"].unique().item()

                    transition_intensities = transition_df_peptide[(transition_df_peptide["transition_id"]==transition) & (transition_df_peptide["data"]==1)]["array"].item()
                    transition_rt = transition_df_peptide[(transition_df_peptide["transition_id"]==transition) & (transition_df_peptide["data"]==2)]["array"].item()

                    decoded_transition_intensities = decode_intensity(transition_intensities)
                    decoded_transition_rt = decode_retention_time(transition_rt)

                    l1, = plt.plot(decoded_transition_rt,decoded_transition_intensities, color=clrs[i+1], linestyle="-", alpha=1)
                    plot_lines_1.append(l1)
                    feature_label_1.append(str(transition_annotation))
                    i+=1
                else:
                    break 
            
            legend_1 = plt.legend(plot_lines_1, feature_label_1, title=f"{peptide}", bbox_to_anchor=(1.05,1), loc="upper left",fontsize=7)

    plt.xlim(min(RT_1, RT_2)-30, max(RT_1, RT_2+30))        
    plt.gca().add_artist(legend_1)
    plt.gca().add_artist(legend_2)
    plt.xlabel("Retention time (s)", fontsize=14)
    plt.ylabel("Intensity", fontsize=14)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    #plt.title(f"{peptide_1} and {peptide_2}", loc="left",fontsize=10)
    plt.axhline(y=0, color = "black", linewidth=1)
    plt.tight_layout()
    plt.savefig(f"{feature_1}_{feature_2}_sqmass_raw_unique_phospho_fragments_precursors_chromatogram.png", dpi=300, bbox_inches="tight")   
    plt.clf()    
    plt.close('all')
