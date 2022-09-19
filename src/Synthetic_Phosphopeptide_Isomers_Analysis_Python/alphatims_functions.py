import alphatims.utils
import alphatims.bruker
import alphatims.plotting
import numpy as np
import holoviews as hv
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.signal import find_peaks
from scipy.signal import peak_widths
from scipy.signal import savgol_filter
from scipy.signal import peak_prominences
from export_format_functions import *
import sqlite3

def get_feature_indices_overlay(tims_data, mz, IM_1, IM_2, RT_1, RT_2, mz_tolerance=50, rt_tolerance=15, im_tolerance=0.05):
    mz_slice = slice(mz/(1+mz_tolerance/10**6), mz*(1+mz_tolerance/10**6))
    rt_slice = slice(float(min(RT_1, RT_2) - rt_tolerance), float(max(RT_1, RT_2)) + rt_tolerance)
    im_slice = slice(float(min(IM_1, IM_2) - im_tolerance), float(max(IM_1, IM_2) + im_tolerance)) 
    feature_indices = tims_data[rt_slice, im_slice, 0, mz_slice, "raw"]
    return feature_indices

def get_transition_indices_overlay(tims_data, transition_mz, IM_1, IM_2, RT_1, RT_2, precursor_mz, mz_tolerance=50, rt_tolerance=15, im_tolerance=0.05):
    transition_mz_slice = slice(transition_mz/(1+mz_tolerance/10**6), transition_mz*(1+mz_tolerance/10**6))
    precursor_mz_slice = slice(precursor_mz/(1+mz_tolerance/10**6), precursor_mz*(1+mz_tolerance/10**6))
    rt_slice = slice(float(min(RT_1, RT_2) - rt_tolerance), float(max(RT_1, RT_2)) + rt_tolerance)
    im_slice = slice(float(min(IM_1, IM_2) - im_tolerance), float(max(IM_1, IM_2) + im_tolerance)) 
    transition_indices = tims_data[rt_slice, im_slice, precursor_mz_slice, transition_mz_slice, "raw"]

    return transition_indices

def get_feature_indices(tims_data, mz, IM, RT,  mz_tolerance=50, rt_tolerance=15, im_tolerance=0.05):
    mz_slice = slice(mz/(1+mz_tolerance/10**6), mz*(1+mz_tolerance/10**6))
    rt_slice = slice(RT - rt_tolerance, RT + rt_tolerance)
    im_slice = slice(IM - im_tolerance, IM + im_tolerance) 

    feature_indices = tims_data[rt_slice, im_slice, 0, mz_slice, "raw"]

    return feature_indices

def precursor_mobilogram(feature_im_values, feature_intensities, smoothed_intensity, fullpeptidename, charge, smoothed_im_peak, sequence, feature, output_folder):
     plt.style.use("seaborn")  
     plt.plot(feature_im_values, feature_intensities, label="Raw values")
     plt.plot(feature_im_values, smoothed_intensity, label="Smoothed Values") 
     plt.xlabel("Ion Mobility Vs $cm^{-2}$")
     plt.ylabel("Intensity")
     plt.title(f"{fullpeptidename}\nCharge state={charge}\nIM={smoothed_im_peak}", loc="left")
     plt.figlegend()
     plt.savefig(f"{output_folder}/{feature}_{sequence}.png", dpi=300)
     plt.clf()
     plt.close('all')

def get_transition_indices(tims_data, transition_mz, IM, RT, precursor_mz, mz_tolerance=50, rt_tolerance=15, im_tolerance=0.05):
    transition_mz_slice = slice(transition_mz/(1+mz_tolerance/10**6), transition_mz*(1+mz_tolerance/10**6))
    precursor_mz_slice = slice(precursor_mz/(1+mz_tolerance/10**6), precursor_mz*(1+mz_tolerance/10**6))
    rt_slice = slice(RT - rt_tolerance, RT + rt_tolerance)
    im_slice = slice(IM - im_tolerance, IM + im_tolerance) 

    transition_indices = tims_data[rt_slice, im_slice, precursor_mz_slice, transition_mz_slice, "raw"]

    return transition_indices


def export_identification_transition_df(osw_file):
    osw_file = sqlite3.connect(osw_file)
    osw_export_query = '''
                   SELECT
                   FEATURE.ID as feature_id,
                        FEATURE.EXP_RT as feature_retention_time,
                        FEATURE.EXP_IM as feature_IM, 
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
 
 
def calculate_mass(mz, annotation):
    charge = float(annotation[-1])
    mz = float(mz)
    mass = charge*mz
    return mass

#transition_df.loc[:, "transition_mass"] = np.vectorize(calculate_mass)(transition_df["transition_mz"], transition_df["transition_annotation"])
