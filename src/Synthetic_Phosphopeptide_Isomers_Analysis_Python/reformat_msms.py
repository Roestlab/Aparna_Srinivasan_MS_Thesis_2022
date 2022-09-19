import pandas as pd
import numpy as np 
import click 
from spectral_lib_functions import *

@click.command()
@click.option('--msms-file', help="full path to msms.txt file", required=True, type=click.Path(exists=True))
@click.option('--msms-output', help="full path to msms output file", required=True, type=str)
def main(msms_file, msms_output):
    """Reformat the phosphorylation modification format in the msms file for diapysef compatibility."""
    msms = pd.read_csv(msms_file, sep = "\t", dtype = {"Modified sequence" : str})

    msms["Modified sequence"] = msms["Modified sequence"].str.replace("Phospho \\(STY\\)", "ph")
    
    msms.to_csv(msms_output, sep="\t", index = False)

if __name__ == "__main__":
    main() 