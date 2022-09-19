import re

def modification_mq_unimod(modifiedsequence):
    """Reformat the modified peptide sequence column in the MQ evidence.txt output for the probability extraction function"""
    modifiedsequence = modifiedsequence.replace("(Phospho (STY))", "(UniMod:21)")
    modifiedsequence = modifiedsequence.replace("(Oxidation (M))", "(UniMod:35)")
    modifiedsequence = modifiedsequence.replace("(Acetyl (Protein N-term))", "(UniMod:1)")
    modifiedsequence = modifiedsequence.replace("_", "")

    return modifiedsequence

#these two functions are the same?
def modification_only_phospho(modifiedsequence):
    """ only keep phosphorylation modification from modified sequence in unimod format"""
    modifiedsequence = modifiedsequence.replace("(UniMod:35)", "")
    modifiedsequence = modifiedsequence.replace("(UniMod:4)",  "")
    modifiedsequence = modifiedsequence.replace("(UniMod:1)",  "")
    return modifiedsequence

def prob_sequence_formatting(modifiedsequence):
    modifiedsequence = modifiedsequence.replace("(Phospho (STY))", "p")
    modifiedsequence = modifiedsequence.replace("(Oxidation (M))", "")
    modifiedsequence = modifiedsequence.replace("(Acetyl (Protein N-term))", "")
    modifiedsequence = modifiedsequence.replace("_", "")

    return modifiedsequence

def phosphorylation_probabilities(probability_sequence):
    """Extract phosphorylation probabilities from the Phospho mods column, then pad the list so it is the same length as the residues list """
    
    padded_seq = []
    
    for i in range(len(probability_sequence)):
        padded_seq.append(probability_sequence[i])
        if i == len(probability_sequence) -1:
            break
        else:
            if probability_sequence[i] == "S" or probability_sequence[i] == "T" or probability_sequence[i] == "Y":
                if probability_sequence[i+1] == "(":
                    continue
                else:
                    padded_seq.append("(0.0)")

    padded_seq = ''.join(padded_seq)
    probabilities_list = [float(x) for x in re.findall("(\d\\.{0,1}\d{0,})", padded_seq)]
       
    return(probabilities_list)

def extract_probability(prob_threshold, mod_sequence, prob_sequence, n_phospho):
    """Determine whether each phosphorylated residue (as reported in the modified sequence) has a probability passing the user defined threshold"""
    if n_phospho==0:
        return 0
    else: 
        phos_residues = re.findall("[STY]p{0,1}", mod_sequence)
        phos_probabilities = phosphorylation_probabilities(prob_sequence)
        phos_position = [i for i in range(len(phos_residues)) if re.search("p", phos_residues[i])]
        phos_probabilities =  [phos_probabilities[i] for i in phos_position] 
        t= 0
        for i in phos_probabilities: 
            if i >= prob_threshold: 
                t+=1
        return t
