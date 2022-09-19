
from alphatims_functions import *
from mobilogram_functions import *
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

            phosphorylated_fragments_1 = fragments_with_phosphorylation(peptide_1)
            phosphorylated_fragments_2 = fragments_with_phosphorylation(peptide_2)

            peptide_1_df = feature_df.loc[feature_df["feature_id"]==feature_1]
            peptide_1_df = peptide_1_df.loc[peptide_1_df["transition_annotation_nocharge"].isin(phosphorylated_fragments_1)]
            peptide_1_phosphorylated_fragments = len(peptide_1_df["transition_annotation_nocharge"].unique())

            peptide_2_df = feature_df.loc[feature_df["feature_id"]==feature_2]
            peptide_2_df = peptide_2_df.loc[peptide_2_df["transition_annotation_nocharge"].isin(phosphorylated_fragments_2)]
            peptide_2_phosphorylated_fragments = len(peptide_2_df["transition_annotation_nocharge"].unique())
        

            peptides_df = pd.concat([peptide_1_df, peptide_2_df]) 
            
            peptides_df["n_features"] = (peptides_df.groupby(["transition_annotation_nocharge"])["feature_id"].transform("nunique")) 
            peptides_df = peptides_df.loc[peptides_df["n_features"] != 2] 
            
            peptide_1_unique_phosphorylated_fragments = len(peptides_df.loc[peptides_df["feature_id"]==feature_1]["transition_annotation_nocharge"].unique())
            peptide_2_unique_phosphorylated_fragments = len(peptides_df.loc[peptides_df["feature_id"]==feature_2]["transition_annotation_nocharge"].unique())


            isomer_df.loc[i, "phosphorylated_fragments.x"] = peptide_1_phosphorylated_fragments
            isomer_df.loc[i, "phosphorylated_fragments.y"] = peptide_2_phosphorylated_fragments
            isomer_df.loc[i, "unique_phosphorylated_fragments.x"] = peptide_1_unique_phosphorylated_fragments
            isomer_df.loc[i, "unique_phosphorylated_fragments.y"] = peptide_2_unique_phosphorylated_fragments

            mobilogram_overlay_precursors(feature_1, feature_2,feature_df, tims_data)
            mobilogram_overlay_precursors_unique_phosphorylated_transitions_annotated(feature_1, feature_2, feature_df, tims_data, feature_1_pep, feature_2_pep, retention_time_delta)
            
        isomer_df.to_csv(f"{raw_data_filename}_isomer_df.tsv", sep="\t", index=False)

if __name__ == "__main__":
    main()

