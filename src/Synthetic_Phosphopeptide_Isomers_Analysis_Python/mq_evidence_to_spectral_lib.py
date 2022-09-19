import pandas as pd
import numpy as np 
import click 
from spectral_lib_functions import *
from pathlib import Path

@click.command()
@click.option('--evidence-file', help="full path to evidence file", required=True, type=click.Path(exists=True))
@click.option('--output-file', default="evidence_output.txt", type=str, required = True, help="Full path and name of filtered evidence file. Will be written to provided path.")
@click.option("--min-probability", default=0.75, required = True, help= "probability cutoff for phosphorylation sites", type=float)
@click.option("--ground-truth", help="ground truth table", type=click.Path(exists=True), required = True)
def main(evidence_file, output_file, min_probability, ground_truth):
    """Read in evidence.txt file from MQ output, annotate phosphorylation probabilities, reformat the modifications and write out the filtered evidence file with iRT peptides included 
    that can be used for diapysef library generation. Also writes out an annotated evidence file with Site localisation information and annotations of whether peptide was found
    in the correct sample. """

    data_folder = Path(output_file)
    data_folder = data_folder.parent

    # read in files
    evidence_table = pd.read_csv(evidence_file, sep = "\t", dtype = {"Modified sequence" : str, "Phospho (STY) Probabilities": str, "Phospho (STY)": int, "id":int})

    ground_truth_table = pd.read_csv(ground_truth, sep = "\t")

    probability_threshold = min_probability
    
    #Remove contaminants and reverse hits from evidence file
    evidence_table = evidence_table.loc[(evidence_table["Potential contaminant"] != "+") & (evidence_table["Reverse"] != "+")]
    evidence_table = evidence_table.dropna(subset=["Proteins"])     
    evidence_table["FullPeptideName_phospho"] = evidence_table["Modified sequence"].apply(lambda x: modification_only_phospho(modification_mq_unimod(x)))

    
    # keep irt peptides separate to be appended to filtered table later 
    irt_table = evidence_table.loc[evidence_table["Proteins"].str.contains("iRT_protein")]

    # reformat ground truth modifications to be consistent with modification format for this step 
    ground_truth_table = ground_truth_table.rename(columns = {"U_pool": "Experiment"})   # rename for easy merging with probability_Table
    ground_truth_table["ground_truth_upool"] = "True"
    ground_truth_table["ground_truth_sequence_upool"] = "True"

    #annotate evidence table 
    evidence_table = evidence_table.merge(ground_truth_table[["FullPeptideName_phospho","Experiment", "ground_truth_upool"]], how="left", on = ["FullPeptideName_phospho", "Experiment"]) #all peptides that are found in the correct upool will have ground_truth annotated as True, while the ones that are not will be annotated as NA
    evidence_table["ground_truth_upool"] = evidence_table["ground_truth_upool"].fillna("False") # all modified peptides that are not found in correct pool are "False" 
    evidence_table = evidence_table.merge(ground_truth_table[["Sequence","Experiment", "ground_truth_sequence_upool"]], how="left", on = ["Sequence", "Experiment"]) #all peptides that are found in the correct upool will have ground_truth annotated as True, while the ones that are not will be annotated as NA
    evidence_table["ground_truth_sequence_upool"] = evidence_table["ground_truth_sequence_upool"].fillna("False") # all peptide sequence that are not found in correct pool are "False" 
    evidence_table["ground_truth_sequence"] = evidence_table["Sequence"].apply(lambda x: "True" if x in ground_truth_table["Sequence"].values else "False")
    evidence_table["ground_truth_localisation"] = evidence_table["FullPeptideName_phospho"].apply(lambda x: "True" if x in ground_truth_table["FullPeptideName_phospho"].values else "False")
    
    # table with unmodified peptides 
    unmodified_peptides_evidence_table = evidence_table.loc[(evidence_table["Phospho (STY)"] == 0) & (evidence_table["ground_truth_sequence_upool"] == "True")] 
    unmodified_peptides_evidence_table = unmodified_peptides_evidence_table.drop(columns=["ground_truth_upool", "ground_truth_sequence_upool", "ground_truth_sequence", "ground_truth_localisation", "FullPeptideName_phospho"],axis=1 )

    #create separate table to do phosphorylation site probability extraction 
    # this table only contains modified peptides which is why the unmodified peptides were left out 
    probability_table = evidence_table[["id", "Modified sequence", "Phospho (STY) Probabilities", "Phospho (STY)", "Experiment", "Sequence"]]  # during maxquant run, experiment was set as the upool sample ID
    probability_table = probability_table.dropna(subset=["Modified sequence", "Phospho (STY) Probabilities", "Phospho (STY)"], how ="any")
    probability_table["Modified sequence"] = probability_table["Modified sequence"].apply(lambda x: prob_sequence_formatting(x))

    #calculate how many phosphorylation sites on the peptide pass threshold or not 
    probability_table.loc[:,"SiteProb"] = np.vectorize(extract_probability)(probability_threshold, probability_table["Modified sequence"], probability_table["Phospho (STY) Probabilities"], probability_table["Phospho (STY)"]) 

    #evidence table for downstream analysis 
    evidence_annotated = pd.merge(evidence_table, probability_table[["id", "SiteProb"]], on="id", how="left") 
    evidence_annotated["SiteProb"] = evidence_annotated["SiteProb"].fillna(0)
    evidence_annotated.to_csv(f"{data_folder}/evidence_annotated.tsv", sep="\t", index=False)  # write the evidence file annotated with Site probability information, and whether modified peptide is in correct u pool etc. 

    # merge probability information with evidence table, format table
    # add back iRT Peptides and unmodified sequences

    evidence_probability_info = evidence_annotated.loc[(evidence_annotated["SiteProb"] == evidence_annotated["Phospho (STY)"]) & (evidence_annotated["ground_truth_upool"] == "True") ]  #this includes the unmodified peptides 
    evidence_probability_info["Modified sequence"] = evidence_probability_info["Modified sequence"].str.replace("Phospho \\(STY\\)", "ph", regex=True)
    
    evidence_final = pd.concat([irt_table, evidence_probability_info, unmodified_peptides_evidence_table])
    evidence_final = evidence_final.rename(columns = {"1/K0" : "IonMobilityIndexK0"})
    evidence_final = evidence_final.drop(columns=["SiteProb", "ground_truth_upool", "ground_truth_sequence_upool", "ground_truth_sequence", "ground_truth_localisation", "FullPeptideName_phospho"], axis=1)

    evidence_final.to_csv(output_file, sep="\t", index = False)   # this is the final evidence file to be used for transition list creation 


if __name__ == "__main__":
    main()  