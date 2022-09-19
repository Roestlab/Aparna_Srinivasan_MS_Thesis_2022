import pandas as pd
import numpy as np 
import click 
from spectral_lib_functions import *

@click.command()
@click.option('--evidence-file', help="full path to evidence file", required=True, type=click.Path(exists=True))
@click.option('--output-file', default = "evidence_default.tsv", required=True, type=str, help="Full path and name of output file. Will be written to provided path")
@click.argument("columns", nargs=-1)
def main(evidence_file, output_file,columns):
    """ Format input evidence file for downstream data analysis """

    # read in table and remove contaminants and reverse hits
    evidence_table = pd.read_csv(evidence_file, sep = "\t", dtype = {"Modified sequence" : str})
    
    # add optional columns
    req_columns = ["Sequence", "Modified sequence", "Proteins", "Experiment", "Charge", "m/z", "Calibrated retention time", "Calibrated retention time start", "Calibrated retention time finish",
    "1/K0", "1/K0 length", "PEP", "Intensity"] 
    
    #create total columns list
    optional_columns = list(columns)

    for column in optional_columns:
        req_columns.append(column)

    #subset columns
    evidence_table = evidence_table[req_columns]

    # add method and accumulation time

    evidence_table["method"] = "DDA_pasef"
    evidence_table["ramp_time"] = 100

    # rename columns
    evidence_table = evidence_table.rename(columns = {"Modified sequence": "FullPeptideName", "Experiment": "Pool", "m/z": "mz", "Calibrated retention time": "RT", "Calibrated retention time start": "leftWidth",
    "Calibrated retention time finish": "rightWidth", "1/K0": "IM", "1/K0 length": "IM_length", "Proteins": "Protein"})

    #format modifications 
    evidence_table["FullPeptideName"] = evidence_table["FullPeptideName"].apply(lambda x: modification_mq_unimod(x))
    evidence_table["FullPeptideName_phospho"] = evidence_table["FullPeptideName"].apply(lambda x: modification_only_phospho(x))

    # write final table
    evidence_table.to_csv(output_file, sep="\t", index=False)

if __name__== "__main__":
    main()