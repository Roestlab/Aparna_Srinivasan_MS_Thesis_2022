from alphatims_functions import * 
from raw_im_extraction_function import *
import click 
import os

@click.command()
@click.option('--path-to-raw-data', help="Path to raw MS data",  type=click.Path(exists=True))
@click.option('--output-folder', help="Output folder for IM plots and dataframe", type=str)
@click.option('--feature-dataframe', help="Dataframe containing features to be extracted", type=click.Path(exists=True))
def main(path_to_raw_data, output_folder, feature_dataframe):
    raw_data_filename = str(os.path.basename(path_to_raw_data))
    raw_data_filename = raw_data_filename[:-2]
    raw_data = alphatims.bruker.TimsTOF(path_to_raw_data)

    feature_df = pd.read_csv(feature_dataframe, sep="\t")
    feature_df = feature_df.loc[feature_df["filename"].str.contains(raw_data_filename) == True] 

    im_extraction_df = raw_ion_mobility_data_extraction(raw_data, feature_df, output_folder)

    im_extraction_df.to_csv(f"{output_folder}/{raw_data_filename}_alphatims_IM_extraction.tsv", sep="\t", index=False)

if __name__ == "__main__":
    main()