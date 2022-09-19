from alphatims_functions import *

def mobilogram_precursor_transitions(feature, feature_df, tims_data):
    
    mobilogram_df = feature_df.loc[feature_df["feature_id"] == feature]
    mobilogram_df = mobilogram_df.loc[mobilogram_df.groupby("transition_annotation")["transition_identifying"].idxmax()]

    im_array = tims_data.mobility_values
    
    peptide = mobilogram_df[mobilogram_df["feature_id"]==feature]["FullPeptideName"].unique().item()
    precursor_mz = mobilogram_df[mobilogram_df["feature_id"]==feature]["precursor_mz"].unique().item()
    rt = mobilogram_df[mobilogram_df["feature_id"]==feature]["feature_retention_time"].unique().item()
    IM = mobilogram_df[mobilogram_df["feature_id"]==feature]["feature_IM"].unique().item()
            
    feature_indices = get_feature_indices(tims_data, precursor_mz, IM, rt)
    feature_intensities = tims_data.bin_intensities(feature_indices, ["mobility_values"])

    clrs_1 = sns.color_palette('Blues', n_colors=len(set(mobilogram_df.loc[mobilogram_df["transition_detecting"] ==1]["transition_id"])))
    clrs_2 = sns.color_palette('Greens', n_colors=len(set(mobilogram_df.loc[mobilogram_df["transition_identifying"] ==1]["transition_id"])))

    i=0
    j=0
    
    plt.figure(figsize=(10,8))
    plt.style.use("seaborn")

    non_zeros = np.flatnonzero(feature_intensities)

    if len(non_zeros) > 0:
        start = max(0, non_zeros[0] - 1)
        end = non_zeros[-1] + 2
        im_values = im_array[start: end]
        feature_intensities = feature_intensities[start: end]

        smoothed_intensity = savgol_filter(feature_intensities, 9,2)
        plt.plot(im_values, smoothed_intensity, color='black', linestyle="-", alpha=1, label=f"{peptide}")

    for i in range(len(mobilogram_df)): 
            transition = mobilogram_df.iloc[i]["transition_id"].item()
            precursor_mz = mobilogram_df[mobilogram_df["transition_id"]==transition]["precursor_mz"].item()
            transition_mz = mobilogram_df[mobilogram_df["transition_id"]==transition]["transition_mz"].item()
            rt = mobilogram_df[mobilogram_df["transition_id"]==transition]["feature_retention_time"].item()
            IM = mobilogram_df[mobilogram_df["transition_id"]==transition]["feature_IM"].item()
            transition_annotation = mobilogram_df[mobilogram_df["transition_id"]==transition]["transition_annotation"].item()
            
            transition_indices = get_transition_indices(tims_data, transition_mz, IM, rt, precursor_mz)
            transition_intensities = tims_data.bin_intensities(transition_indices, ["mobility_values"])

            non_zeros = np.flatnonzero(transition_intensities)
            if len(non_zeros) <= 10:
                continue

            else:
                start = max(0, non_zeros[0] - 1)
                end = non_zeros[-1] + 2
                im_values = im_array[start: end]
                transition_intensities = transition_intensities[start: end]

                smoothed_intensity = savgol_filter(transition_intensities, 9,2)
    
                if mobilogram_df.loc[mobilogram_df["transition_id"] ==transition]["transition_detecting"].item() ==1:
                    plt.plot(im_values, smoothed_intensity, color=clrs_1[i], linestyle="--", alpha=1, label=f"{transition_annotation}")
                    i+=1
                else:
                    plt.plot(im_values, smoothed_intensity, color=clrs_2[i], linestyle="--", alpha=1, label=f"{transition_annotation}")
                    j+=1

    plt.legend(bbox_to_anchor=(1.05,1), loc="upper left", fontsize=7)
    plt.xlabel("Ion Mobility Vs $cm^{-2}$")
    plt.ylabel("Intensity")
    plt.axhline(y=0, color = "black", linewidth=1)
    plt.tight_layout()
    plt.title(f"{peptide}", loc="left", fontsize=10)
    plt.savefig(f"{feature}_mobilogram_precursor_fragments.png", dpi=300, bbox_inches="tight")   
    plt.clf()
    plt.close('all')

def mobilogram_precursor_phosphorylated_transitions(feature, feature_df, tims_data):
 
    mobilogram_df = feature_df.loc[feature_df["feature_id"] == feature]
    peptide = mobilogram_df["FullPeptideName"].unique().item()   #this could be wrong
    mobilogram_df = mobilogram_df.loc[mobilogram_df.groupby("transition_annotation")["transition_identifying"].idxmax()]
    mobilogram_df = mobilogram_df.loc[mobilogram_df["transition_annotation_nocharge"].isin(fragments_with_phosphorylation(peptide))]

    im_array = tims_data.mobility_values
    
    precursor_mz = mobilogram_df[mobilogram_df["feature_id"]==feature]["precursor_mz"].unique().item()

    rt = mobilogram_df[mobilogram_df["feature_id"]==feature]["feature_retention_time"].unique().item()
    IM = mobilogram_df[mobilogram_df["feature_id"]==feature]["feature_IM"].unique().item()
            
    feature_indices = get_feature_indices(tims_data, precursor_mz, IM, rt)
    feature_intensities = tims_data.bin_intensities(feature_indices, ["mobility_values"])

    clrs = sns.color_palette('Greens', n_colors=len(mobilogram_df))

    i=0
    
    plt.figure(figsize=(10,8))
    plt.style.use("seaborn")

    non_zeros = np.flatnonzero(feature_intensities)

    if len(non_zeros) > 0:
        start = max(0, non_zeros[0] - 1)
        end = non_zeros[-1] + 2
        im_values = im_array[start: end]
        feature_intensities = feature_intensities[start: end]

        smoothed_intensity = savgol_filter(feature_intensities, 9,2)
        plt.plot(im_values, smoothed_intensity, color='black', linestyle="-", alpha=1, label=f"{peptide}")

    for i in range(len(mobilogram_df)): 
            transition = mobilogram_df.iloc[i]["transition_id"].item()
            precursor_mz = mobilogram_df[mobilogram_df["transition_id"]==transition]["precursor_mz"].item()
            transition_mz = mobilogram_df[mobilogram_df["transition_id"]==transition]["transition_mz"].item()
            rt = mobilogram_df[mobilogram_df["transition_id"]==transition]["feature_retention_time"].item()
            IM = mobilogram_df[mobilogram_df["transition_id"]==transition]["feature_IM"].item()
            transition_annotation = mobilogram_df[mobilogram_df["transition_id"]==transition]["transition_annotation"].item()
            
            transition_indices = get_transition_indices(tims_data, transition_mz, IM, rt, precursor_mz)
            transition_intensities = tims_data.bin_intensities(transition_indices, ["mobility_values"])

            non_zeros = np.flatnonzero(transition_intensities)
            if len(non_zeros) <= 10:
                continue

            else:
                start = max(0, non_zeros[0] - 1)
                end = non_zeros[-1] + 2
                im_values = im_array[start: end]
                transition_intensities = transition_intensities[start: end]

                smoothed_intensity = savgol_filter(transition_intensities, 9,2)
    
                plt.plot(im_values, smoothed_intensity, color=clrs[i], linestyle="--", alpha=1, label=f"{transition_annotation}")
        

    plt.legend(bbox_to_anchor=(1.05,1), loc="upper left", fontsize=7)
    plt.xlabel("Ion Mobility Vs $cm^{-2}$")
    plt.ylabel("Intensity")
    plt.axhline(y=0, color = "black", linewidth=1)
    plt.tight_layout()
    plt.title(f"{peptide}", loc="left", fontsize=10)
    plt.savefig(f"{feature}_mobilogram_precursor_phosphorylated_fragments.png", dpi=300, bbox_inches="tight")   
    plt.clf()
    plt.close('all')

def mobilogram_transitions(feature, feature_df, tims_data):
    mobilogram_df = feature_df.loc[feature_df["feature_id"] == feature]
    mobilogram_df = mobilogram_df.loc[mobilogram_df.groupby("transition_annotation")["transition_identifying"].idxmax()]

    im_array = tims_data.mobility_values
    
    peptide = mobilogram_df[mobilogram_df["feature_id"]==feature]["FullPeptideName"].unique().item()
    precursor_mz = mobilogram_df[mobilogram_df["feature_id"]==feature]["precursor_mz"].unique().item()
    rt = mobilogram_df[mobilogram_df["feature_id"]==feature]["feature_retention_time"].unique().item()
    IM = mobilogram_df[mobilogram_df["feature_id"]==feature]["feature_IM"].unique().item()

    clrs_1 = sns.color_palette('Blues', n_colors=len(set(mobilogram_df.loc[mobilogram_df["transition_detecting"] ==1]["transition_id"])))
    clrs_2 = sns.color_palette('Greens', n_colors=len(set(mobilogram_df.loc[mobilogram_df["transition_identifying"] ==1]["transition_id"])))

    i=0
    j=0
    
    plt.figure(figsize=(10,8))
    plt.style.use("seaborn")

    for i in range(len(mobilogram_df)): 
            transition = mobilogram_df.iloc[i]["transition_id"].item()
            transition_mz = mobilogram_df[mobilogram_df["transition_id"]==transition]["transition_mz"].item()
            transition_annotation = mobilogram_df[mobilogram_df["transition_id"]==transition]["transition_annotation"].item()
            
            transition_indices = get_transition_indices(tims_data, transition_mz, IM, rt, precursor_mz)
            transition_intensities = tims_data.bin_intensities(transition_indices, ["mobility_values"])

            non_zeros = np.flatnonzero(transition_intensities)
            if len(non_zeros) <= 10:
                continue

            else:
                start = max(0, non_zeros[0] - 1)
                end = non_zeros[-1] + 2
                im_values = im_array[start: end]
                transition_intensities = transition_intensities[start: end]

                smoothed_intensity = savgol_filter(transition_intensities, 9,2)
    
                if mobilogram_df.loc[mobilogram_df["transition_id"] ==transition]["transition_detecting"].item() ==1:
                    plt.plot(im_values, smoothed_intensity, color=clrs_1[i], linestyle="--", alpha=1, label=f"{transition_annotation}")
                    i+=1
                else:
                    plt.plot(im_values, smoothed_intensity, color=clrs_2[i], linestyle="--", alpha=1, label=f"{transition_annotation}")
                    j+=1

    plt.legend(bbox_to_anchor=(1.05,1), loc="upper left", fontsize=7)
    plt.xlabel("Ion Mobility Vs $cm^{-2}$")
    plt.ylabel("Intensity")
    plt.axhline(y=0, color = "black", linewidth=1)
    plt.tight_layout()
    plt.title(f"{peptide}", loc="left", fontsize=10)
    plt.savefig(f"{feature}_mobilogram_fragments.png", dpi=300, bbox_inches="tight")   
    plt.clf() 
    plt.close('all')

def mobilogram_phosphorylated_transitions(feature, feature_df, tims_data):
    mobilogram_df = feature_df.loc[feature_df["feature_id"] == feature]
    mobilogram_df = mobilogram_df.loc[mobilogram_df.groupby("transition_annotation")["transition_identifying"].idxmax()]
    peptide = mobilogram_df["FullPeptideName"].unique().item() 
    mobilogram_df = mobilogram_df.loc[mobilogram_df["transition_annotation_nocharge"].isin(fragments_with_phosphorylation(peptide))]

    im_array = tims_data.mobility_values
    
      #this could be wrong
    precursor_mz = mobilogram_df[mobilogram_df["feature_id"]==feature]["precursor_mz"].unique().item()
    rt = mobilogram_df[mobilogram_df["feature_id"]==feature]["feature_retention_time"].unique().item()
    IM = mobilogram_df[mobilogram_df["feature_id"]==feature]["feature_IM"].unique().item()
            
    clrs = sns.color_palette('Greens', n_colors=len(mobilogram_df))

    i=0
    
    plt.figure(figsize=(10,8))
    plt.style.use("seaborn")

    for i in range(len(mobilogram_df)): 
            transition = mobilogram_df.iloc[i]["transition_id"].item()
            precursor_mz = mobilogram_df[mobilogram_df["transition_id"]==transition]["precursor_mz"].item()
            transition_mz = mobilogram_df[mobilogram_df["transition_id"]==transition]["transition_mz"].item()
            rt = mobilogram_df[mobilogram_df["transition_id"]==transition]["feature_retention_time"].item()
            IM = mobilogram_df[mobilogram_df["transition_id"]==transition]["feature_IM"].item()
            transition_annotation = mobilogram_df[mobilogram_df["transition_id"]==transition]["transition_annotation"].item()
            
            transition_indices = get_transition_indices(tims_data, transition_mz, IM, rt, precursor_mz)
            transition_intensities = tims_data.bin_intensities(transition_indices, ["mobility_values"])

            non_zeros = np.flatnonzero(transition_intensities)
            if len(non_zeros) <= 10:
                continue

            else:
                start = max(0, non_zeros[0] - 1)
                end = non_zeros[-1] + 2
                im_values = im_array[start: end]
                transition_intensities = transition_intensities[start: end]

                smoothed_intensity = savgol_filter(transition_intensities, 9,2)
    
                plt.plot(im_values, smoothed_intensity, color=clrs[i], linestyle="--", alpha=1, label=f"{transition_annotation}")
        

    plt.legend(bbox_to_anchor=(1.05,1), loc="upper left", fontsize=7)
    plt.xlabel("Ion Mobility Vs $cm^{-2}$")
    plt.ylabel("Intensity")
    plt.axhline(y=0, color = "black", linewidth=1)
    plt.tight_layout()
    plt.title(f"{peptide}", loc="left", fontsize=10)
    plt.savefig(f"{feature}_mobilogram_phosphorylated_fragments.png", dpi=300, bbox_inches="tight")   
    plt.clf()
    plt.close('all')

def mobilogram_overlay_precursors(feature_1, feature_2, feature_df, tims_data):

    mobilogram_df = feature_df.loc[feature_df["feature_id"].isin([feature_1, feature_2])]

    peptides = mobilogram_df["FullPeptideName"].unique()

    im_array = tims_data.mobility_values

    plt.figure(figsize=(10,8))
    for feature in [feature_1, feature_2]:
            peptide = mobilogram_df[mobilogram_df["feature_id"]==feature]["FullPeptideName"].unique().item()
            precursor_mz = mobilogram_df[mobilogram_df["feature_id"]==feature]["precursor_mz"].unique().item()
            rt = mobilogram_df[mobilogram_df["feature_id"]==feature]["feature_retention_time"].unique().item()
            IM = mobilogram_df[mobilogram_df["feature_id"]==feature]["library_IM"].unique().item()
            
            feature_indices = get_feature_indices(tims_data, precursor_mz, IM, rt)
            feature_intensities = tims_data.bin_intensities(feature_indices, ["mobility_values"])

            non_zeros = np.flatnonzero(feature_intensities)
            if len(non_zeros) == 0:
                continue
            else:
                start = max(0, non_zeros[0] - 1)
                end = non_zeros[-1] + 2
                im_values = im_array[start: end]
                feature_intensities = feature_intensities[start: end]

                smoothed_intensity = savgol_filter(feature_intensities, 9,2)
           
                if feature == feature_1: 
                    plt.plot(im_values, smoothed_intensity, color="forestgreen", linestyle="-", alpha=1, label=f"{peptide}")
                    plt.text(IM, np.amax(smoothed_intensity) + 100, f"Library IM = {IM}", horizontalalignment='center')
                else: 
                    plt.plot(im_values, smoothed_intensity, color="deepskyblue", linestyle="-", alpha=1, label=f"{peptide}")
                    plt.text(IM, np.amax(smoothed_intensity) + 100, f"Library IM = {IM}", horizontalalignment='center')

        
    plt.style.use("seaborn")
    plt.legend(bbox_to_anchor=(1.05,1), loc="upper left", fontsize=7)
    plt.xlabel("Ion Mobility Vs $cm^{-2}$")
    plt.ylabel("Intensity")
    plt.title(f"{peptides[0]} and {peptides[1]}", loc="left", fontsize=10)
    plt.axhline(y=0, color = "black", linewidth=1)
    plt.tight_layout()
    plt.savefig(f"{feature_1}_{feature_2}_precursors_mobilogram.png", dpi=300)   
    plt.clf()  
    plt.close('all')
    plt.close('all')  

def mobilogram_overlay_precursors_unique_transitions(feature_1, feature_2, feature_df, tims_data):
    
    mobilogram_df = feature_df.loc[feature_df["feature_id"].isin([feature_1, feature_2])]
    mobilogram_df["n_features"] = (mobilogram_df.groupby(["transition_mz", "transition_annotation"])["feature_id"].transform("nunique"))
    mobilogram_df = mobilogram_df.loc[mobilogram_df["n_features"] != 2]
    mobilogram_df = mobilogram_df.loc[mobilogram_df.groupby(["feature_id", "transition_annotation"])["transition_identifying"].idxmax()]

    clrs_1 =  sns.color_palette('Blues', n_colors= len(mobilogram_df.loc[mobilogram_df["feature_id"] == feature_1])+1)
    clrs_2 =  sns.color_palette('Greens', n_colors= len(mobilogram_df.loc[mobilogram_df["feature_id"] == feature_2])+1)

    peptides = mobilogram_df["FullPeptideName"].unique() # order of peptides might not correspond to feature_1 and feature_2 

    im_array = tims_data.mobility_values

    plot_lines_1 = []
    plot_lines_2 = []
    feature_label_1 = []
    feature_label_2 = []

    plt.figure(figsize=(10,8))


    feature_list = np.array([feature_1, feature_2])
    phosphorylation_positions = np.array([phosphorylation_position(peptide_1), phosphorylation_position(peptide_2)])
    phosphorylation_positions = np.argsort(phosphorylation_positions)
    features_sorted = feature_list[phosphorylation_positions]

    for feature in features_sorted:
            peptide = mobilogram_df[mobilogram_df["feature_id"]==feature]["FullPeptideName"].unique().item()
            precursor_mz = mobilogram_df[mobilogram_df["feature_id"]==feature]["precursor_mz"].unique().item()
            rt = mobilogram_df[mobilogram_df["feature_id"]==feature]["feature_retention_time"].unique().item()
            IM = mobilogram_df[mobilogram_df["feature_id"]==feature]["feature_IM"].unique().item()
            
            feature_indices = get_feature_indices(tims_data, precursor_mz, IM, rt)
            feature_intensities = tims_data.bin_intensities(feature_indices, ["mobility_values"])

            non_zeros = np.flatnonzero(feature_intensities)
            if len(non_zeros) == 0:
                continue
            else:
                start = max(0, non_zeros[0] - 1)
                end = non_zeros[-1] + 2
                im_values = im_array[start: end]
                feature_intensities = feature_intensities[start: end]

                smoothed_intensity = savgol_filter(feature_intensities, 9,2)
           
                if feature == feature_1: 
                    l1, = plt.plot(im_values, smoothed_intensity, color="deepskyblue", linestyle="-", alpha=1, label=f"{peptide}_precursor")
                    plot_lines_1.append(l1)
                    feature_label_1.append("precursor")
                else: 
                    l2, = plt.plot(im_values, smoothed_intensity, color="forestgreen", linestyle="-", alpha=1, label=f"{peptide}_precursor")
                    plot_lines_1.append(l2)
                    feature_label_2.append("precursor")
            
            mobilogram_df_peptide = mobilogram_df.loc[mobilogram_df["feature_id"] == feature] 
 
            for i in range(len(mobilogram_df_peptide)):
                
                if len(mobilogram_df_peptide) > 0:

                    transition = mobilogram_df_peptide.iloc[i]["transition_id"].item()
                    transition_mz = mobilogram_df[mobilogram_df["transition_id"]==transition]["transition_mz"].item()
                    transition_annotation = mobilogram_df_peptide[mobilogram_df_peptide["transition_id"]==transition]["transition_annotation"].item()
                        
                    transition_indices = get_transition_indices(tims_data, transition_mz, IM, rt, precursor_mz)
                    transition_intensities = tims_data.bin_intensities(transition_indices, ["mobility_values"])

                    non_zeros = np.flatnonzero(transition_intensities)
                    if len(non_zeros) <= 10:
                        continue

                    else:
                        start = max(0, non_zeros[0] - 1)
                        end = non_zeros[-1] + 2
                        im_values = im_array[start: end]
                        transition_intensities = transition_intensities[start: end]

                        smoothed_intensity = savgol_filter(transition_intensities, 9,2)
                    
                        if feature == feature_1: 
                            l1, = plt.plot(im_values, smoothed_intensity, color=clrs_1[i+1], linestyle="--", alpha=1)
                            plot_lines_1.append(l1)
                            feature_label_1.append(str(transition_annotation))
                        else: 
                            l2, = plt.plot(im_values, smoothed_intensity, color=clrs_2[i+1], linestyle="--", alpha=1, label=f"{transition_annotation}_{peptide}")
                            plot_lines_2.append(l2)
                            feature_label_2.append(str(transition_annotation))
                else: 
                    break

    plt.style.use("seaborn")
    legend_1 = plt.legend(plot_lines_1, feature_label_1, title=f"{peptides[0]}", bbox_to_anchor=(1.05,1), loc="upper left",fontsize=7)
    legend_2 = plt.legend(plot_lines_2, feature_label_2, title=f"{peptides[1]}", bbox_to_anchor=(1.05,0.5), loc="center left",fontsize=7)
    plt.gca().add_artist(legend_1)
    plt.gca().add_artist(legend_2)
    plt.xlabel("Ion Mobility Vs $cm^{-2}$")
    plt.ylabel("Intensity")
    plt.title(f"{peptides[0]} and {peptides[1]}", loc="left", fontsize=10)
    plt.axhline(y=0, color = "black", linewidth=1)
    plt.tight_layout()
    plt.savefig(f"{feature_1}_{feature_2}_unique_fragments_precursor_mobilogram.png", dpi=300, bbox_inches="tight")   
    plt.clf()
    plt.close('all')

def mobilogram_overlay_unique_transitions(feature_1, feature_2, feature_df, tims_data):
    mobilogram_df = feature_df.loc[feature_df["feature_id"].isin([feature_1, feature_2])]
    mobilogram_df["n_features"] = (mobilogram_df.groupby(["transition_mz", "transition_annotation"])["feature_id"].transform("nunique"))
    mobilogram_df = mobilogram_df.loc[mobilogram_df["n_features"] != 2]
    mobilogram_df = mobilogram_df.loc[mobilogram_df.groupby(["feature_id", "transition_annotation"])["transition_identifying"].idxmax()]

    clrs_1 =  sns.color_palette('Blues', n_colors= len(mobilogram_df.loc[mobilogram_df["feature_id"] == feature_1])+1)
    clrs_2 =  sns.color_palette('Greens', n_colors= len(mobilogram_df.loc[mobilogram_df["feature_id"] == feature_2])+1)

    peptides = mobilogram_df["FullPeptideName"].unique()

    im_array = tims_data.mobility_values

    plot_lines_1 = []
    plot_lines_2 = []
    feature_label_1 = []
    feature_label_2 = []

    plt.figure(figsize=(10,8))


    feature_list = np.array([feature_1, feature_2])
    phosphorylation_positions = np.array([phosphorylation_position(peptide_1), phosphorylation_position(peptide_2)])
    phosphorylation_positions = np.argsort(phosphorylation_positions)
    features_sorted = feature_list[phosphorylation_positions]

    for feature in features_sorted:
            
            precursor_mz = mobilogram_df[mobilogram_df["feature_id"]==feature]["precursor_mz"].unique().item()
            rt = mobilogram_df[mobilogram_df["feature_id"]==feature]["feature_retention_time"].unique().item()
            IM = mobilogram_df[mobilogram_df["feature_id"]==feature]["feature_IM"].unique().item()
            
            mobilogram_df_peptide = mobilogram_df.loc[mobilogram_df["feature_id"] == feature] 
 
            for i in range(len(mobilogram_df_peptide)):
                
                if len(mobilogram_df_peptide) > 0:

                    transition = mobilogram_df_peptide.iloc[i]["transition_id"].item()
                    peptide = mobilogram_df_peptide[mobilogram_df_peptide["transition_id"]==transition]["FullPeptideName"].item()
                    transition_mz = mobilogram_df_peptide[mobilogram_df_peptide["transition_id"]==transition]["transition_mz"].item()
                    transition_annotation = mobilogram_df_peptide[mobilogram_df_peptide["transition_id"]==transition]["transition_annotation"].item()
                        
                    transition_indices = get_transition_indices(tims_data, transition_mz, IM, rt, precursor_mz)
                    transition_intensities = tims_data.bin_intensities(transition_indices, ["mobility_values"])

                    non_zeros = np.flatnonzero(transition_intensities)
                    if len(non_zeros) <= 10:
                        continue

                    else:
                        start = max(0, non_zeros[0] - 1)
                        end = non_zeros[-1] + 2
                        im_values = im_array[start: end]
                        transition_intensities = transition_intensities[start: end]

                        smoothed_intensity = savgol_filter(transition_intensities, 9,2)
                    
                        if feature == feature_1: 
                            l1, = plt.plot(im_values, smoothed_intensity, color=clrs_1[i+1], linestyle="--", alpha=1)
                            plot_lines_1.append(l1)
                            feature_label_1.append(str(transition_annotation))
                        else: 
                            l2, = plt.plot(im_values, smoothed_intensity, color=clrs_2[i+1], linestyle="--", alpha=1, label=f"{transition_annotation}_{peptide}")
                            plot_lines_2.append(l2)
                            feature_label_2.append(str(transition_annotation))
                else: 
                    break

    plt.style.use("seaborn")
    legend_1 = plt.legend(plot_lines_1, feature_label_1, title=f"{peptides[0]}", bbox_to_anchor=(1.05,1), loc="upper left",fontsize=7)
    legend_2 = plt.legend(plot_lines_2, feature_label_2, title=f"{peptides[1]}", bbox_to_anchor=(1.05,0.5), loc="center left",fontsize=7)
    plt.gca().add_artist(legend_1)
    plt.gca().add_artist(legend_2)
    plt.xlabel("Ion Mobility Vs $cm^{-2}$")
    plt.ylabel("Intensity")
    plt.title(f"{peptides[0]} and {peptides[1]}", loc="left", fontsize=10)
    plt.axhline(y=0, color = "black", linewidth=1)
    plt.tight_layout()
    plt.savefig(f"{feature_1}_{feature_2}_unique_fragments_mobilogram.png", dpi=300, bbox_inches="tight")   
    plt.clf()
    plt.close('all')


def mobilogram_overlay_precursors_unique_phosphorylated_transitions(feature_1, feature_2, feature_df, tims_data):

    peptide_df_1 = feature_df.loc[feature_df["feature_id"] == feature_1]
    peptide_1 = feature_df.loc[feature_df["feature_id"] == feature_1]["FullPeptideName"].unique().item()
    peptide_df_1 = peptide_df_1.loc[peptide_df_1["transition_annotation_nocharge"].isin(fragments_with_phosphorylation(peptide_1))]

    peptide_df_2 = feature_df.loc[feature_df["feature_id"] == feature_2]
    peptide_2 = feature_df.loc[feature_df["feature_id"] == feature_2]["FullPeptideName"].unique().item()
    peptide_df_2 = peptide_df_2.loc[peptide_df_2["transition_annotation_nocharge"].isin(fragments_with_phosphorylation(peptide_2))]

    peptide_df = pd.concat([peptide_df_1, peptide_df_2], axis=0)


    peptide_df["n_features"] = (peptide_df.groupby(["transition_annotation_nocharge"])["feature_id"].transform("nunique")) 
    peptide_df = peptide_df.loc[peptide_df["n_features"] != 2] 
    peptide_df = peptide_df.loc[peptide_df.groupby(["feature_id", "transition_annotation", "transition_mz"])["transition_identifying"].idxmax()]

    clrs_1 =  sns.color_palette('Blues', n_colors= len(peptide_df.loc[peptide_df["feature_id"] == feature_1])+1)
    clrs_2 =  sns.color_palette('Greens', n_colors= len(peptide_df.loc[peptide_df["feature_id"] == feature_2])+1)

    im_array = tims_data.mobility_values

    plt.figure(figsize=(10,8))
    plt.style.use("seaborn")


    plot_lines_1 = []
    plot_lines_2 = []
    feature_label_1 = []
    feature_label_2 = []

    plt.figure(figsize=(10,8))


    feature_list = np.array([feature_1, feature_2])
    phosphorylation_positions = np.array([phosphorylation_position(peptide_1), phosphorylation_position(peptide_2)])
    phosphorylation_positions = np.argsort(phosphorylation_positions)
    features_sorted = feature_list[phosphorylation_positions]


    for feature in features_sorted:
        peptide = feature_df.loc[feature_df["feature_id"]==feature]["FullPeptideName"].unique().item()
        precursor_mz = feature_df.loc[feature_df["feature_id"]==feature]["precursor_mz"].unique().item()
        rt = feature_df.loc[feature_df["feature_id"]==feature]["feature_retention_time"].unique().item()
        IM = feature_df.loc[feature_df["feature_id"]==feature]["library_IM"].unique().item()

        feature_indices = get_feature_indices(tims_data, precursor_mz, IM, rt)
        feature_intensities = tims_data.bin_intensities(feature_indices, ["mobility_values"])

        non_zeros = np.flatnonzero(feature_intensities)
        if len(non_zeros) == 0:
            continue
        else:
            start = max(0, non_zeros[0] - 1)
            end = non_zeros[-1] + 2
            im_values = im_array[start: end]
            feature_intensities = feature_intensities[start: end]

            smoothed_intensity = savgol_filter(feature_intensities, 9,2)
            smoothed_intensity = smoothed_intensity*0.1
            if feature == feature_1: 
                l1, = plt.plot(im_values, smoothed_intensity, color="deepskyblue", linestyle="-", alpha=1)
                plt.text(IM, np.amax(smoothed_intensity) + 100, f"Library IM = {IM}", horizontalalignment='center')
                plot_lines_1.append(l1)
                feature_label_1.append("precursor")

            else: 
                l2, = plt.plot(im_values, smoothed_intensity, color="forestgreen", linestyle="-", alpha=1, label=f"{peptide}_precursor")
                plt.text(IM, np.amax(smoothed_intensity) + 100, f"Library IM = {IM}", horizontalalignment='center')
                plot_lines_2.append(l2)
                feature_label_2.append("precursor")

        transition_df_peptide = peptide_df.loc[peptide_df["feature_id"] == feature] 


        for i in range(len(transition_df_peptide)):
            if len(transition_df_peptide) > 0:
                transition = transition_df_peptide.iloc[i]["transition_id"]
                transition_mz = transition_df_peptide[transition_df_peptide["transition_id"]==transition]["transition_mz"].item()
                transition_annotation = transition_df_peptide[transition_df_peptide["transition_id"]==transition]["transition_annotation"].item()

                transition_indices = get_transition_indices(tims_data, transition_mz, IM, rt, precursor_mz)
                transition_intensities = tims_data.bin_intensities(transition_indices, ["mobility_values"])

                non_zeros = np.flatnonzero(transition_intensities)
                if len(non_zeros) <= 10:
                    continue

                else:
                    start = max(0, non_zeros[0] - 1)
                    end = non_zeros[-1] + 2
                    im_values = im_array[start: end]
                    transition_intensities = transition_intensities[start: end]

                    smoothed_intensity = savgol_filter(transition_intensities, 9,2)

                    if feature == feature_1: 
                        l1, = plt.plot(im_values, smoothed_intensity, color=clrs_1[i+1], linestyle="--", alpha=1)
                        plot_lines_1.append(l1)
                        feature_label_1.append(str(transition_annotation))
                    else: 
                        l2, = plt.plot(im_values, smoothed_intensity, color=clrs_2[i+1], linestyle="--", alpha=1)
                        plot_lines_2.append(l2)
                        feature_label_2.append(str(transition_annotation))
            else:
                break 

    plt.style.use("seaborn")
    legend_1 = plt.legend(plot_lines_1, feature_label_1, title=f"{peptide_1}", bbox_to_anchor=(1.05,1), loc="upper left",fontsize=7)
    legend_2 = plt.legend(plot_lines_2, feature_label_2, title=f"{peptide_2}", bbox_to_anchor=(1.05,0.5), loc="center left",fontsize=7)
    plt.gca().add_artist(legend_1)
    plt.gca().add_artist(legend_2)
    plt.xlabel("Ion Mobility Vs $cm^{-2}$")
    plt.ylabel("Intensity")
    plt.title(f"{peptide_1} and {peptide_2}", loc="left",fontsize=10)
    plt.axhline(y=0, color = "black", linewidth=1)
    plt.tight_layout()
    plt.savefig(f"{feature_1}_{feature_2}_unique_phospho_fragments_precursors_mobilogram.png", dpi=300, bbox_inches="tight")   
    plt.clf()    
    plt.close('all')

def mobilogram_overlay_unique_phosphorylated_transition(feature_1, feature_2, feature_df, tims_data):

    peptide_df_1 = feature_df.loc[feature_df["feature_id"] == feature_1]
    peptide_1 = feature_df.loc[feature_df["feature_id"] == feature_1]["FullPeptideName"].unique().item()
    peptide_df_1 = peptide_df_1.loc[peptide_df_1["transition_annotation_nocharge"].isin(fragments_with_phosphorylation(peptide_1))]

    peptide_df_2 = feature_df.loc[feature_df["feature_id"] == feature_2]
    peptide_2 = feature_df.loc[feature_df["feature_id"] == feature_2]["FullPeptideName"].unique().item()
    peptide_df_2 = peptide_df_2.loc[peptide_df_2["transition_annotation_nocharge"].isin(fragments_with_phosphorylation(peptide_2))]

    peptide_df = pd.concat([peptide_df_1, peptide_df_2], axis=0)


    peptide_df["n_features"] = (peptide_df.groupby(["transition_annotation"])["feature_id"].transform("nunique"))
    peptide_df = peptide_df.loc[peptide_df["n_features"] != 2] #this doesn't get rid of the case that they might have the same transition but at different charge states the m/z is different. Maybe just go with annotation for now? 
    peptide_df = peptide_df.loc[peptide_df.groupby(["feature_id", "transition_annotation", "transition_mz"])["transition_identifying"].idxmax()]

    clrs_1 =  sns.color_palette('Blues', n_colors= len(peptide_df.loc[peptide_df["feature_id"] == feature_1])+1)
    clrs_2 =  sns.color_palette('Greens', n_colors= len(peptide_df.loc[peptide_df["feature_id"] == feature_2])+1)

    im_array = tims_data.mobility_values

    plt.figure(figsize=(10,8))
    plt.style.use("seaborn")


    plot_lines_1 = []
    plot_lines_2 = []
    feature_label_1 = []
    feature_label_2 = []

    plt.figure(figsize=(10,8))


    feature_list = np.array([feature_1, feature_2])
    phosphorylation_positions = np.array([phosphorylation_position(peptide_1), phosphorylation_position(peptide_2)])
    phosphorylation_positions = np.argsort(phosphorylation_positions)
    features_sorted = feature_list[phosphorylation_positions]



    for feature in features_sorted:
        peptide = feature_df.loc[feature_df["feature_id"]==feature]["FullPeptideName"].unique().item()
        precursor_mz = feature_df.loc[feature_df["feature_id"]==feature]["precursor_mz"].unique().item()
        rt = feature_df.loc[feature_df["feature_id"]==feature]["feature_retention_time"].unique().item()
        IM = feature_df.loc[feature_df["feature_id"]==feature]["feature_IM"].unique().item()

        transition_df_peptide = peptide_df.loc[peptide_df["feature_id"] == feature] 

        for i in range(len(transition_df_peptide)):
            if len(transition_df_peptide) > 0:
                transition = transition_df_peptide.iloc[i]["transition_id"]
                transition_mz = transition_df_peptide[transition_df_peptide["transition_id"]==transition]["transition_mz"].item()
                transition_annotation = transition_df_peptide[transition_df_peptide["transition_id"]==transition]["transition_annotation"].item()

                transition_indices = get_transition_indices(tims_data, transition_mz, IM, rt, precursor_mz)
                transition_intensities = tims_data.bin_intensities(transition_indices, ["mobility_values"])

                non_zeros = np.flatnonzero(transition_intensities)
                if len(non_zeros) <= 10:
                    continue

                else:
                    start = max(0, non_zeros[0] - 1)
                    end = non_zeros[-1] + 2
                    im_values = im_array[start: end]
                    transition_intensities = transition_intensities[start: end]

                    smoothed_intensity = savgol_filter(transition_intensities, 9,2)

                    if feature == feature_1: 
                        l1, = plt.plot(im_values, smoothed_intensity, color=clrs_1[i+1], linestyle="--", alpha=1, label=f"{transition_annotation}_{peptide}")
                        plot_lines_1.append(l1)
                        feature_label_1.append(str(transition_annotation))
                    else: 
                        l2, = plt.plot(im_values, smoothed_intensity, color=clrs_2[i+1], linestyle="--", alpha=1, label=f"{transition_annotation}_{peptide}")
                        plot_lines_2.append(l2)
                        feature_label_2.append(str(transition_annotation))
            else:
                break 

    plt.style.use("seaborn")
    legend_1 = plt.legend(plot_lines_1, feature_label_1, title=f"{peptide_1}", bbox_to_anchor=(1.05,1), loc="upper left",fontsize=7)
    legend_2 = plt.legend(plot_lines_2, feature_label_2, title=f"{peptide_2}", bbox_to_anchor=(1.05,0.5), loc="center left",fontsize=7)
    plt.gca().add_artist(legend_1)
    plt.gca().add_artist(legend_2)
    plt.xlabel("Ion Mobility Vs $cm^{-2}$")
    plt.ylabel("Intensity")
    plt.title(f"{peptide_1} and {peptide_2}", loc="left", fontsize=10)
    plt.axhline(y=0, color = "black", linewidth=1)
    plt.tight_layout()
    plt.savefig(f"{feature_1}_{feature_2}_unique_phospho_fragments_mobilogram.png", dpi=300,  bbox_inches="tight")   
    plt.clf()
    plt.close('all')    

def mobilogram_overlay_all_phosphorylated_transition(feature_1, feature_2, feature_df, tims_data):

    peptide_df_1 = feature_df.loc[feature_df["feature_id"] == feature_1]
    peptide_1 = feature_df.loc[feature_df["feature_id"] == feature_1]["FullPeptideName"].unique().item()
    peptide_df_1 = peptide_df_1.loc[peptide_df_1["transition_annotation_nocharge"].isin(fragments_with_phosphorylation(peptide_1))]

    peptide_df_2 = feature_df.loc[feature_df["feature_id"] == feature_2]
    peptide_2 = feature_df.loc[feature_df["feature_id"] == feature_2]["FullPeptideName"].unique().item()
    peptide_df_2 = peptide_df_2.loc[peptide_df_2["transition_annotation_nocharge"].isin(fragments_with_phosphorylation(peptide_2))]

    peptide_df = pd.concat([peptide_df_1, peptide_df_2], axis=0)

    peptide_df = peptide_df.loc[peptide_df.groupby(["feature_id", "transition_annotation", "transition_mz"])["transition_identifying"].idxmax()]

    clrs_1 =  sns.color_palette('Blues', n_colors= len(peptide_df.loc[peptide_df["feature_id"] == feature_1])+1)
    clrs_2 =  sns.color_palette('Greens', n_colors= len(peptide_df.loc[peptide_df["feature_id"] == feature_2])+1)

    im_array = tims_data.mobility_values

    plt.figure(figsize=(10,8))
    plt.style.use("seaborn")


    plot_lines_1 = []
    plot_lines_2 = []
    feature_label_1 = []
    feature_label_2 = []

    plt.figure(figsize=(10,8))


    feature_list = np.array([feature_1, feature_2])
    phosphorylation_positions = np.array([phosphorylation_position(peptide_1), phosphorylation_position(peptide_2)])
    phosphorylation_positions = np.argsort(phosphorylation_positions)
    features_sorted = feature_list[phosphorylation_positions]
    
    for feature in features_sorted:
        peptide = feature_df.loc[feature_df["feature_id"]==feature]["FullPeptideName"].unique().item()
        precursor_mz = feature_df.loc[feature_df["feature_id"]==feature]["precursor_mz"].unique().item()
        rt = feature_df.loc[feature_df["feature_id"]==feature]["feature_retention_time"].unique().item()
        IM = feature_df.loc[feature_df["feature_id"]==feature]["feature_IM"].unique().item()

        transition_df_peptide = peptide_df.loc[peptide_df["feature_id"] == feature] 

        for i in range(len(transition_df_peptide)):
            if len(transition_df_peptide) > 0:
                transition = transition_df_peptide.iloc[i]["transition_id"]
                transition_mz = transition_df_peptide[transition_df_peptide["transition_id"]==transition]["transition_mz"].item()
                transition_annotation = transition_df_peptide[transition_df_peptide["transition_id"]==transition]["transition_annotation"].item()

                transition_indices = get_transition_indices(tims_data, transition_mz, IM, rt, precursor_mz)
                transition_intensities = tims_data.bin_intensities(transition_indices, ["mobility_values"])

                non_zeros = np.flatnonzero(transition_intensities)
                if len(non_zeros) <= 10:
                    continue

                else:
                    start = max(0, non_zeros[0] - 1)
                    end = non_zeros[-1] + 2
                    im_values = im_array[start: end]
                    transition_intensities = transition_intensities[start: end]

                    smoothed_intensity = savgol_filter(transition_intensities, 9,2)

                    if feature == feature_1: 
                        l1, = plt.plot(im_values, smoothed_intensity, color=clrs_1[i+1], linestyle="--", alpha=1, label=f"{transition_annotation}_{peptide}")
                        plot_lines_1.append(l1)
                        feature_label_1.append(str(transition_annotation))
                    else: 
                        l2, = plt.plot(im_values, smoothed_intensity, color=clrs_2[i+1], linestyle="--", alpha=1, label=f"{transition_annotation}_{peptide}")
                        plot_lines_2.append(l2)
                        feature_label_2.append(str(transition_annotation))
            else:
                break 

    plt.style.use("seaborn")
    legend_1 = plt.legend(plot_lines_1, feature_label_1, title=f"{peptide_1}", bbox_to_anchor=(1.05,1), loc="upper left",fontsize=7)
    legend_2 = plt.legend(plot_lines_2, feature_label_2, title=f"{peptide_2}", bbox_to_anchor=(1.05,0.5), loc="center left",fontsize=7)
    plt.gca().add_artist(legend_1)
    plt.gca().add_artist(legend_2)
    plt.xlabel("Ion Mobility Vs $cm^{-2}$")
    plt.ylabel("Intensity")
    plt.title(f"{peptide_1} and {peptide_2}", loc="left", fontsize=10)
    plt.axhline(y=0, color = "black", linewidth=1)
    plt.tight_layout()
    plt.savefig(f"{feature_1}_{feature_2}_phospho_fragments_mobilogram.png", dpi=300,  bbox_inches="tight")   
    plt.clf()
    plt.close('all')    

def mobilogram_overlay_precursors_unique_phosphorylated_transitions_annotated(feature_1, feature_2, feature_df, tims_data, feature_1_pep, feature_2_pep, retention_time_delta):
    peptide_df_1 = feature_df.loc[feature_df["feature_id"] == feature_1]
    peptide_1 = feature_df.loc[feature_df["feature_id"] == feature_1]["FullPeptideName"].unique().item()

    peptide_df_2 = feature_df.loc[feature_df["feature_id"] == feature_2]
    peptide_2 = feature_df.loc[feature_df["feature_id"] == feature_2]["FullPeptideName"].unique().item()
    
    peptide_1_fragments = fragments_with_phosphorylation(peptide_1)
    peptide_2_fragments = fragments_with_phosphorylation(peptide_2)
    peptide_1_unique_fragments = set(peptide_1_fragments).difference(peptide_2_fragments)
    peptide_2_unique_fragments = set(peptide_2_fragments).difference(peptide_1_fragments)

    peptide_df_2 = peptide_df_2.loc[peptide_df_2["transition_annotation_nocharge"].isin(peptide_2_unique_fragments)]
    peptide_df_1 = peptide_df_1.loc[peptide_df_1["transition_annotation_nocharge"].isin(peptide_1_unique_fragments)]


    peptide_df = pd.concat([peptide_df_1, peptide_df_2], axis=0)

    peptide_df = peptide_df.loc[peptide_df.groupby(["feature_id", "transition_annotation", "transition_mz"])["transition_identifying"].idxmax()]

    clrs_1 =  sns.color_palette('Blues', n_colors= len(peptide_df.loc[peptide_df["feature_id"] == feature_1])+1)
    clrs_2 =  sns.color_palette('Greens', n_colors= len(peptide_df.loc[peptide_df["feature_id"] == feature_2])+1)

    im_array = tims_data.mobility_values

    plt.figure(figsize=(10,8))
    plt.style.use("seaborn")


    plot_lines_1 = []
    plot_lines_2 = []
    feature_label_1 = []
    feature_label_2 = []


    feature_list = np.array([feature_1, feature_2])
    phosphorylation_positions = np.array([phosphorylation_position(peptide_1), phosphorylation_position(peptide_2)])
    phosphorylation_positions = np.argsort(phosphorylation_positions)
    features_sorted = feature_list[phosphorylation_positions]

    for feature in features_sorted:
        peptide = feature_df.loc[feature_df["feature_id"]==feature]["FullPeptideName"].unique().item()
        precursor_mz = feature_df.loc[feature_df["feature_id"]==feature]["precursor_mz"].unique().item()
        rt = feature_df.loc[feature_df["feature_id"]==feature]["feature_retention_time"].unique().item()
        IM = feature_df.loc[feature_df["feature_id"]==feature]["library_IM"].unique().item()

        feature_indices = get_feature_indices(tims_data, precursor_mz, IM, rt)
        feature_intensities = tims_data.bin_intensities(feature_indices, ["mobility_values"])

        non_zeros = np.flatnonzero(feature_intensities)
        if len(non_zeros) == 0:
            continue
        else:
            start = max(0, non_zeros[0] - 1)
            end = non_zeros[-1] + 2
            im_values = im_array[start: end]
            feature_intensities = feature_intensities[start: end]

            smoothed_intensity = savgol_filter(feature_intensities, 9,2)
            smoothed_intensity = smoothed_intensity*0.1
            if feature == feature_1: 
                l1, = plt.plot(im_values, smoothed_intensity, color="deepskyblue", linestyle="-", alpha=1)
                plt.axvline(x=IM, color="deepskyblue", linewidth=1, linestyle=":")
                #plt.text(IM, np.amax(smoothed_intensity), f"Library IM = {IM}", horizontalalignment='left', color="deepskyblue")
                plot_lines_1.append(l1)
                feature_label_1.append("precursor")

            else: 
                l2, = plt.plot(im_values, smoothed_intensity, color="forestgreen", linestyle="-", alpha=1, label=f"{peptide}_precursor")
                plt.axvline(x=IM, color="forestgreen", linewidth=1, linestyle=":")
                #plt.text(IM, np.amax(smoothed_intensity), f"Library IM = {IM}", horizontalalignment='right', color="forestgreen")
                plot_lines_2.append(l2)
                feature_label_2.append("precursor")

        transition_df_peptide = peptide_df.loc[peptide_df["feature_id"] == feature] 


        for i in range(len(transition_df_peptide)):
            if len(transition_df_peptide) > 0:
                transition = transition_df_peptide.iloc[i]["transition_id"]
                transition_mz = transition_df_peptide[transition_df_peptide["transition_id"]==transition]["transition_mz"].item()
                transition_annotation = transition_df_peptide[transition_df_peptide["transition_id"]==transition]["transition_annotation"].item()

                transition_indices = get_transition_indices(tims_data, transition_mz, IM, rt, precursor_mz)
                transition_intensities = tims_data.bin_intensities(transition_indices, ["mobility_values"])

                non_zeros = np.flatnonzero(transition_intensities)
                if len(non_zeros) <= 10:
                    continue

                else:
                    start = max(0, non_zeros[0] - 1)
                    end = non_zeros[-1] + 2
                    im_values = im_array[start: end]
                    transition_intensities = transition_intensities[start: end]

                    smoothed_intensity = savgol_filter(transition_intensities, 9,2)

                    if feature == feature_1: 
                        l1, = plt.plot(im_values, smoothed_intensity, color=clrs_1[i+1], linestyle="--", alpha=1)
                        plot_lines_1.append(l1)
                        feature_label_1.append(str(transition_annotation))
                    else: 
                        l2, = plt.plot(im_values, smoothed_intensity, color=clrs_2[i+1], linestyle="--", alpha=1)
                        plot_lines_2.append(l2)
                        feature_label_2.append(str(transition_annotation))
            else:
                break 

    plt.style.use("seaborn")
    legend_1 = plt.legend(plot_lines_1, feature_label_1, title=f"{peptide_1}", bbox_to_anchor=(1.05,1), loc="upper left",fontsize=7)
    legend_2 = plt.legend(plot_lines_2, feature_label_2, title=f"{peptide_2}", bbox_to_anchor=(1.05,0.5), loc="center left",fontsize=7)
    plt.gca().add_artist(legend_1)
    plt.gca().add_artist(legend_2)
    plt.xlabel("Ion Mobility Vs $cm^{-2}$")
    plt.ylabel("Intensity")
    plt.title(f"Retention time difference = {retention_time_delta}s \n IPF PEP scores: {peptide_1} = {feature_1_pep}, {peptide_2} = {feature_2_pep}", loc="left",fontsize=10)
    plt.axhline(y=0, color = "black", linewidth=1)
    plt.tight_layout()
    plt.savefig(f"{feature_1}_{feature_2}_unique_phospho_fragments_precursors_mobilogram.png", dpi=300, bbox_inches="tight")   
    plt.clf()    
    plt.close('all')

def mobilogram_overlay_precursors_unique_phosphorylated_transitions_im_recalibration_annotated(feature_1, feature_2, feature_df, tims_data, feature_1_pep, feature_2_pep, retention_time_delta, transition_im_x, transition_im_y):
    peptide_df_1 = feature_df.loc[feature_df["feature_id"] == feature_1]
    peptide_1 = feature_df.loc[feature_df["feature_id"] == feature_1]["FullPeptideName"].unique().item()

    peptide_df_2 = feature_df.loc[feature_df["feature_id"] == feature_2]
    peptide_2 = feature_df.loc[feature_df["feature_id"] == feature_2]["FullPeptideName"].unique().item()
    
    peptide_1_fragments = fragments_with_phosphorylation(peptide_1)
    peptide_2_fragments = fragments_with_phosphorylation(peptide_2)
    peptide_1_unique_fragments = set(peptide_1_fragments).difference(peptide_2_fragments)
    peptide_2_unique_fragments = set(peptide_2_fragments).difference(peptide_1_fragments)

    peptide_df_2 = peptide_df_2.loc[peptide_df_2["transition_annotation_nocharge"].isin(peptide_2_unique_fragments)]
    peptide_df_1 = peptide_df_1.loc[peptide_df_1["transition_annotation_nocharge"].isin(peptide_1_unique_fragments)]


    peptide_df = pd.concat([peptide_df_1, peptide_df_2], axis=0)

    peptide_df = peptide_df.loc[peptide_df.groupby(["feature_id", "transition_annotation", "transition_mz"])["transition_identifying"].idxmax()]

    clrs_1 =  sns.color_palette('Blues', n_colors= len(peptide_df.loc[peptide_df["feature_id"] == feature_1])+1)
    clrs_2 =  sns.color_palette('Greens', n_colors= len(peptide_df.loc[peptide_df["feature_id"] == feature_2])+1)

    im_array = tims_data.mobility_values

    plt.figure(figsize=(10,8))
    plt.style.use("seaborn")


    plot_lines_1 = []
    plot_lines_2 = []
    feature_label_1 = []
    feature_label_2 = []

    
    feature_list = np.array([feature_1, feature_2])
    phosphorylation_positions = np.array([phosphorylation_position(peptide_1), phosphorylation_position(peptide_2)])
    phosphorylation_positions = np.argsort(phosphorylation_positions)
    features_sorted = feature_list[phosphorylation_positions]

    for feature in features_sorted:
        peptide = feature_df.loc[feature_df["feature_id"]==feature]["FullPeptideName"].unique().item()
        precursor_mz = feature_df.loc[feature_df["feature_id"]==feature]["precursor_mz"].unique().item()
        rt = feature_df.loc[feature_df["feature_id"]==feature]["feature_retention_time"].unique().item()
        IM = feature_df.loc[feature_df["feature_id"]==feature]["library_IM"].unique().item()

        feature_indices = get_feature_indices(tims_data, precursor_mz, IM, rt)
        feature_intensities = tims_data.bin_intensities(feature_indices, ["mobility_values"])

        non_zeros = np.flatnonzero(feature_intensities)
        if len(non_zeros) == 0:
            continue
        else:
            start = max(0, non_zeros[0] - 1)
            end = non_zeros[-1] + 2
            im_values = im_array[start: end]
            feature_intensities = feature_intensities[start: end]

            smoothed_intensity = savgol_filter(feature_intensities, 25,2)
            smoothed_intensity = smoothed_intensity*0.1
            if feature == feature_1: 
                l1, = plt.plot(im_values, smoothed_intensity, color="deepskyblue", linestyle="-", alpha=1)
                if np.isnan(transition_im_x):
                    plt.axvline(x=IM, color="deepskyblue", linewidth=1, linestyle=":")
                else:
                    plt.axvline(x=transition_im_x, color="deepskyblue", linewidth=1, linestyle=":")
                
                plot_lines_1.append(l1)
                feature_label_1.append("precursor")

            else: 
                l2, = plt.plot(im_values, smoothed_intensity, color="forestgreen", linestyle="-", alpha=1, label=f"{peptide}_precursor")

                if np.isnan(transition_im_y):
                    plt.axvline(x=IM, color="forestgreen", linewidth=1, linestyle=":")
                else: 
                    plt.axvline(x=transition_im_y, color="forestgreen", linewidth=1, linestyle=":")
                
                plot_lines_2.append(l2)
                feature_label_2.append("precursor")

        transition_df_peptide = peptide_df.loc[peptide_df["feature_id"] == feature] 


        for i in range(len(transition_df_peptide)):
            if len(transition_df_peptide) > 0:
                transition = transition_df_peptide.iloc[i]["transition_id"]
                transition_mz = transition_df_peptide[transition_df_peptide["transition_id"]==transition]["transition_mz"].item()
                transition_annotation = transition_df_peptide[transition_df_peptide["transition_id"]==transition]["transition_annotation"].item()

                transition_indices = get_transition_indices(tims_data, transition_mz, IM, rt, precursor_mz)
                transition_intensities = tims_data.bin_intensities(transition_indices, ["mobility_values"])

                non_zeros = np.flatnonzero(transition_intensities)
                if len(non_zeros) <= 10:
                    continue

                else:
                    start = max(0, non_zeros[0] - 1)
                    end = non_zeros[-1] + 2
                    im_values = im_array[start: end]
                    transition_intensities = transition_intensities[start: end]

                    smoothed_intensity = savgol_filter(transition_intensities, 25,2)

                    if feature == feature_1: 
                        l1, = plt.plot(im_values, smoothed_intensity, color=clrs_1[i+1], linestyle="--", alpha=1)
                        plot_lines_1.append(l1)
                        feature_label_1.append(str(transition_annotation))
                    else: 
                        l2, = plt.plot(im_values, smoothed_intensity, color=clrs_2[i+1], linestyle="--", alpha=1)
                        plot_lines_2.append(l2)
                        feature_label_2.append(str(transition_annotation))
            else:
                break 

    plt.style.use("seaborn")
    legend_1 = plt.legend(plot_lines_1, feature_label_1, title=f"{peptide_1}", bbox_to_anchor=(1.05,1), loc="upper left",fontsize=7)
    legend_2 = plt.legend(plot_lines_2, feature_label_2, title=f"{peptide_2}", bbox_to_anchor=(1.05,0.5), loc="center left",fontsize=7)
    plt.gca().add_artist(legend_1)
    plt.gca().add_artist(legend_2)
    plt.xlabel("Ion Mobility Vs $cm^{-2}$")
    plt.ylabel("Intensity")
    plt.title(f"Retention time difference = {retention_time_delta}s \n IPF PEP scores: {peptide_1} = {feature_1_pep}, {peptide_2} = {feature_2_pep}", loc="left",fontsize=10)
    plt.axhline(y=0, color = "black", linewidth=1)
    plt.tight_layout()
    plt.savefig(f"{feature_1}_{feature_2}_unique_phospho_fragments_precursors_im_recalibration_mobilogram.png", dpi=300, bbox_inches="tight")   
    plt.clf()    
    plt.close('all')