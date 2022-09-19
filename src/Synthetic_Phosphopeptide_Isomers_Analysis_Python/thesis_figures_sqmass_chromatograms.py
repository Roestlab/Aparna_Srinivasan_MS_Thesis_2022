from sqmass_chromatograms import *
from export_format_functions import *
import click
import re
import sqlite3


@click.command()
@click.option('--path-to-raw-data', help="Path to raw data file",  type=click.Path(exists=True), required=True)
@click.option('--data-folder', help="Path to Scored Pyprophet Folder", required=True)
@click.option('--chromatograms-folder', help="Path to Chromatograms Folder", required=True)
@click.option('--isomer-dataframe', help="Path to isomer dataframe",  type=click.Path(exists=True), required=True)
def main(path_to_raw_data, data_folder, chromatograms_folder, isomer_dataframe):
    raw_data_filename = str(os.path.basename(path_to_raw_data))
    raw_data_filename = raw_data_filename[:-2]

    isomer_df = pd.read_csv(isomer_dataframe, sep="\t")
    isomer_df = isomer_df.loc[isomer_df["filename.x"].str.contains(raw_data_filename)]
    isomer_df = isomer_df.reset_index()
    isomer_df = isomer_df.drop(columns={"index"})

    if len(isomer_df) > 0:
        osw_0 = f"{data_folder}/{raw_data_filename}_0.osw"
        osw_1 = f"{data_folder}/{raw_data_filename}_1.osw"

        sqmass_0 = f"{chromatograms_folder}/{raw_data_filename}_0.sqMass"
        sqmass_1 = f"{chromatograms_folder}/{raw_data_filename}_1.sqMass"

        chromatogram_transitions_1 = export_transition_chromatogram_df(sqmass_0)
        chromatogram_transitions_2 = export_transition_chromatogram_df(sqmass_1)

        chromatogram_precursors_1 = export_precursor_chromatogram_df(sqmass_0)
        chromatogram_precursors_2 = export_precursor_chromatogram_df(sqmass_1)

        feature_df_1 = export_all_transitions(osw_0)
        feature_df_2 = export_all_transitions(osw_1)

        chromatogram_transition_df_1 = merged_transition_df(feature_df_1, chromatogram_transitions_1)
        chromatogram_transition_df_2 = merged_transition_df(feature_df_2, chromatogram_transitions_2)
        chromatogram_transitions = pd.concat([chromatogram_transition_df_1, chromatogram_transition_df_2],axis=0)


        chromatogram_precursor_df_1 = merged_precursor_df(feature_df_1, chromatogram_precursors_1)
        chromatogram_precursor_df_2 = merged_precursor_df(feature_df_2, chromatogram_precursors_2)
        chromatogram_precursors = pd.concat([chromatogram_precursor_df_1, chromatogram_precursor_df_2],axis=0)

        for i in range(len(isomer_df)):
            feature_1 = isomer_df.iloc[i]["feature_id.x"]
            feature_2 = isomer_df.iloc[i]["feature_id.y"]

            chromatogram_overlay_precursors_unique_phosphorylated_transitions_annotated_nosmoothing(feature_1, feature_2, chromatogram_transitions, chromatogram_precursors)



if __name__ == "__main__":
    main()

