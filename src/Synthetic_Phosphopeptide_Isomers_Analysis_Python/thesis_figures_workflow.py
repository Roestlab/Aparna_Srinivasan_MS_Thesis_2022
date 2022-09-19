
# for co-eluting isomer pairs 

from alphatims_functions import *
from mobilogram_functions_nosmoothing import *
from export_format_functions import *
import click
import re
import sqlite3
import os 


@click.command()
@click.option('--path-to-raw-data', help="Path to raw data file",  type=click.Path(exists=True), required=True)
@click.option('--data-folder', help="Path to Scored Pyprophet Folder", required=True)
@click.option('--isomer-dataframe', help="Path to isomer dataframe",  type=click.Path(exists=True), required=True)
def main(path_to_raw_data, data_folder, isomer_dataframe):
    raw_data_filename = str(os.path.basename(path_to_raw_data))
    raw_data_filename = raw_data_filename[:-2]

    isomer_df = pd.read_csv(isomer_dataframe, sep="\t")
    isomer_df = isomer_df.loc[isomer_df["filename.x"].str.contains(raw_data_filename)]
    isomer_df = isomer_df.reset_index()
    isomer_df = isomer_df.drop(columns={"index"})

    if len(isomer_df) > 0:

        tims_data = alphatims.bruker.TimsTOF(path_to_raw_data)

        entries = os.listdir(data_folder)
        osw_list = [entry for entry in entries if entry.endswith(".osw") and bool(re.match(raw_data_filename, entry))]
        df_list = []
        for osw in osw_list:
            df = export_all_transitions(f"{data_folder}/{osw}")
            df_list.append(df) 
        
        feature_df = pd.concat(df_list, ignore_index=True)


        for i in range(len(isomer_df)):
            feature_1 = isomer_df.iloc[i]["feature_id.x"]
            peptide_1 = isomer_df.iloc[i]["FullPeptideName.x"]
            feature_2 = isomer_df.iloc[i]["feature_id.y"]
            peptide_2 = isomer_df.iloc[i]["FullPeptideName.y"]

            feature_1_pep = isomer_df.iloc[i]["ipf_pep.x"]
            feature_2_pep = isomer_df.iloc[i]["ipf_pep.y"]
            retention_time_delta = isomer_df.iloc[i]["delta_RT"]

            mobilogram_overlay_precursors_nosmoothing(feature_1, feature_2,feature_df, tims_data)
            mobilogram_overlay_precursors_unique_phosphorylated_transitions_annotated_nosmoothing(feature_1, feature_2, feature_df, tims_data, feature_1_pep, feature_2_pep, retention_time_delta)
   
            

if __name__ == "__main__":
    main()

