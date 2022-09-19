import pandas as pd
import numpy as np 
import click
from export_format_functions import *
import sqlite3

@click.command()
@click.option('--input-file', help="full path to input file. , will be written to pyprophet output directory", default="pyprophet_formatted_results.tsv", type=click.Path(exists=True))
@click.option('--output-file', help="full path and name of output file", required=True)
def main(input_file, output_file):
    '''Export results from Scored pyprophet osw file. '''
    
    add_mapping_table(input_file) # add unimod codename mapping table to the scored pyprophet file 

    results_df = custom_pyprophet_export(input_file)
    
    results_df = results_df[["feature_id", "Sequence", "FullPeptideName", "precursor_charge", "precursor_mz", "RT", "IM", "Intensity", "leftWidth", "rightWidth", "ms2_pep", "ms2_m_score", "ipf_pep", "ipf_m_score", "peak_group_rank"]]

    results_df = results_df.rename(columns = {"precursor_charge": "Charge", "precursor_mz":"mz"})


    results_df.to_csv(output_file, sep="\t", index = False)

if __name__ == "__main__":
    main() 



