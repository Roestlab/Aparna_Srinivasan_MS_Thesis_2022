import pandas as pd
import numpy as np 
from pathlib import Path
import os
import click
from export_format_functions import *

@click.command()
@click.option('--output-file', help="name of formatted output file, will be written to input directory", required=True)
@click.option('--pool/--no-pool', help = "True or False", default=True)
@click.option('--ramp-time/--no-ramp-time', help="True or false", default=True)
@click.option('--openswath-workflow-number', help="Workflow identification number", type=int)
@click.option('--pyprophet-workflow-number', help="Workflow identification number", type=int)
@click.option('--library', help="spectral library", type=str)
@click.argument('input_files_directory')
def main(pool,ramp_time, output_file, library, openswath_workflow_number, pyprophet_workflow_number, input_files_directory):
    '''    Concatenate data tables into one table and add "meta data" '''
    data_folder = Path(input_files_directory)
    entries = os.listdir(data_folder)

    df_list = []
    for entry in entries: 
        df = pd.read_csv(data_folder/entry, sep = "\t")
        df["filename"] = entry
        df_list.append(df) 
        
    results = pd.concat(df_list, axis=0, ignore_index=True)
    
    results["FullPeptideName_phospho"] = results.apply(lambda row: modification_only_phospho(row["FullPeptideName"]), axis =1)

    results["osw_workflow_id"] = openswath_workflow_number

    results["pyprophet_workflow_id"] = pyprophet_workflow_number
    
    results["library"] = library 

    if pool:
        results.loc[:,"pool"] = np.vectorize(get_pool_name)(results["filename"])

    if ramp_time:
        results.loc[:,"ramp_time"] = np.vectorize(get_ramp_time)(results["filename"])

    results.to_csv(f"{data_folder}/{output_file}", sep="\t", index=False)

if __name__ == "__main__":
    main()