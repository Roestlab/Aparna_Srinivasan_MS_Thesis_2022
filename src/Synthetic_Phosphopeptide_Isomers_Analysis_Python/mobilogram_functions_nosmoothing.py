from alphatims_functions import *
from export_format_functions import * 

def mobilogram_overlay_precursors_unique_phosphorylated_transitions_annotated_nosmoothing(feature_1, feature_2, feature_df, tims_data, feature_1_pep, feature_2_pep, retention_time_delta):

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

    im_array = tims_data.mobility_values

    plt.figure(figsize=(10,8))
    plt.style.use("seaborn")

    plot_lines_1 = []
    plot_lines_2 = []
    feature_label_1 = []
    feature_label_2 = []

    plt.figure(figsize=(10,8))

    RT_1 = peptide_df.loc[peptide_df["feature_id"] == feature_1]["feature_retention_time"].unique().item()
    RT_2 = peptide_df.loc[peptide_df["feature_id"] == feature_2]["feature_retention_time"].unique().item()

    IM_1 = peptide_df.loc[peptide_df["feature_id"] == feature_1]["library_IM"].unique().item()
    IM_2 = peptide_df.loc[peptide_df["feature_id"] == feature_2]["library_IM"].unique().item()

    feature_list = np.array([feature_1, feature_2])
    phosphorylation_positions = np.array([phosphorylation_position(peptide_1), phosphorylation_position(peptide_2)])
    phosphorylation_positions = np.argsort(phosphorylation_positions)
    features_sorted = feature_list[phosphorylation_positions]
    features_sorted = np.flip(features_sorted)

    

    for j in range(2):
        feature = features_sorted[j]
        if j ==1: 
            clrs = sns.color_palette('Blues', n_colors= len(peptide_df.loc[peptide_df["feature_id"] == feature])+1)
        else:
            clrs = sns.color_palette('Greens', n_colors= len(peptide_df.loc[peptide_df["feature_id"] == feature])+1)
    
        precursor_mz = feature_df.loc[feature_df["feature_id"]==feature]["precursor_mz"].unique().item()
    
        transition_df_peptide = peptide_df.loc[peptide_df["feature_id"] == feature] 

        feature_indices = get_feature_indices_overlay(tims_data, precursor_mz, IM_1, IM_2, RT_1, RT_2, im_tolerance=0.1)
        feature_intensities = tims_data.bin_intensities(feature_indices, ["mobility_values"])

        non_zeros = np.flatnonzero(feature_intensities)
        if len(non_zeros) == 0:
            continue
        else:
            start = max(0, non_zeros[0] - 1)
            end = non_zeros[-1] + 2
            im_values = im_array[start: end]
            feature_intensities = feature_intensities[start: end]
            feature_intensities = feature_intensities*0.1

            if feature == feature_1: 
                l1, = plt.plot(im_values, feature_intensities, color="black", linestyle="-", alpha=0.8)
                #plt.text(IM, np.amax(feature_intensities) + 100, f"Library IM = {IM}", horizontalalignment='center')
                plot_lines_1.append(l1)
                feature_label_1.append("precursor")

            else: 
                l2, = plt.plot(im_values, feature_intensities, color="black", linestyle="-", alpha=0.8)
                #plt.text(IM, np.amax(feature_intensities) + 100, f"Library IM = {IM}", horizontalalignment='center')
                plot_lines_2.append(l2)
                feature_label_2.append("precursor")

        for i in range(len(transition_df_peptide)):
            if len(transition_df_peptide) > 0:
                transition = transition_df_peptide.iloc[i]["transition_id"]
                transition_mz = transition_df_peptide[transition_df_peptide["transition_id"]==transition]["transition_mz"].item()
                transition_annotation = transition_df_peptide[transition_df_peptide["transition_id"]==transition]["transition_annotation"].item()

                transition_indices = get_transition_indices_overlay(tims_data, transition_mz, IM_1, IM_2, RT_1, RT_2, precursor_mz, im_tolerance=0.1)
                transition_intensities = tims_data.bin_intensities(transition_indices, ["mobility_values"])

                non_zeros = np.flatnonzero(transition_intensities)
                if len(non_zeros) <= 10:
                    continue

                else:
                    start = max(0, non_zeros[0] - 1)
                    end = non_zeros[-1] + 2
                    im_values = im_array[start: end]
                    transition_intensities = transition_intensities[start: end]

                    if feature == feature_1: 
                        l1, = plt.plot(im_values, transition_intensities, color=clrs[i+1], linestyle="--", alpha=1)
                        plot_lines_1.append(l1)
                        feature_label_1.append(str(transition_annotation))
                    else: 
                        l2, = plt.plot(im_values, transition_intensities, color=clrs[i+1], linestyle="--", alpha=1)
                        plot_lines_2.append(l2)
                        feature_label_2.append(str(transition_annotation))
            else:
                break 
        
    feature_label_1 = feature_label_1[0:21]
    feature_label_2 = feature_label_2[0:21]

    plt.style.use("seaborn")
    legend_1 = plt.legend(plot_lines_1, feature_label_1, title=f"{peptide_1}", bbox_to_anchor=(1.05,1), loc="upper left",fontsize=7)
    legend_2 = plt.legend(plot_lines_2, feature_label_2, title=f"{peptide_2}", bbox_to_anchor=(1.05,0.20), loc="center left",fontsize=7)
    plt.gca().add_artist(legend_1)
    plt.gca().add_artist(legend_2)
    plt.xlabel("Ion Mobility Vs $cm^{-2}$", fontsize=14)
    plt.ylabel("Intensity", fontsize=14)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    #plt.title(f"Retention time difference = {retention_time_delta}s \n IPF PEP scores: {peptide_1} = {feature_1_pep}, {peptide_2} = {feature_2_pep}", loc="left",fontsize=10)
    plt.axhline(y=0, color = "black", linewidth=1)
    plt.tight_layout()
    plt.savefig(f"{feature_1}_{feature_2}_raw_unique_phospho_fragments_precursors_mobilogram.png", dpi=300, bbox_inches="tight")   
    plt.clf()    
    plt.close('all')

def chromatogram_overlay_precursors_nosmoothing(feature_1, feature_2, feature_df, tims_data):

    mobilogram_df = feature_df.loc[feature_df["feature_id"].isin([feature_1, feature_2])]
    rt_array = tims_data.rt_values

    plt.figure(figsize=(10,8))
    plt.style.use("seaborn")
    precursor_mz = mobilogram_df[mobilogram_df["feature_id"]==feature_1]["precursor_mz"].unique().item()
    
    RT_1 = mobilogram_df[mobilogram_df["feature_id"]==feature_1]["feature_retention_time"].unique().item()
    RT_2 = mobilogram_df[mobilogram_df["feature_id"]==feature_2]["feature_retention_time"].unique().item()

    peptide_1 = mobilogram_df[mobilogram_df["feature_id"]==feature_1]["FullPeptideName"].unique().item()
    peptide_2 = mobilogram_df[mobilogram_df["feature_id"]==feature_2]["FullPeptideName"].unique().item()
   
    feature_indices = get_feature_indices_overlay(tims_data, precursor_mz, 1, 1, RT_1, RT_2, rt_tolerance=30, im_tolerance=0.6)   # don't want to consider ion mobility while doing chromatogram extraction
    feature_intensities = tims_data.bin_intensities(feature_indices, ["rt_values"])

    non_zeros = np.flatnonzero(feature_intensities)
    start = max(0, non_zeros[0] - 1)
    end = non_zeros[-1] + 2
    rt_values = rt_array[start: end]
    feature_intensities = feature_intensities[start: end]

    feature_list = np.array([feature_1, feature_2])
    phosphorylation_positions = np.array([phosphorylation_position(peptide_1), phosphorylation_position(peptide_2)])
    phosphorylation_positions = np.argsort(phosphorylation_positions)
    features_sorted = feature_list[phosphorylation_positions]
    features_sorted = np.flip(features_sorted)


    plt.plot(rt_values, feature_intensities, color="black", linestyle="-", alpha=1)      

    for i in range(2):
        feature = features_sorted[i]
        peptide = mobilogram_df[mobilogram_df["feature_id"]==feature]["FullPeptideName"].unique().item()
        RT = mobilogram_df[mobilogram_df["feature_id"]==feature]["feature_retention_time"].unique().item()   
        if i ==1:
            rt_color = "deepskyblue" 
            plt.axvline(x=RT, color=rt_color, linewidth=1, linestyle=":", label=f"{peptide}")
        else:     
            rt_color = "forestgreen"
            plt.axvline(x=RT, color=rt_color, linewidth=1, linestyle=":", label=f"{peptide}")

    plt.legend(bbox_to_anchor=(1.05,1), loc="upper left", fontsize=7)
    plt.xlabel("Retention time (s)")
    plt.ylabel("Intensity")
    plt.title(f"{peptide_1} and {peptide_2}", loc="right", fontsize=10)
    plt.axhline(y=0, color = "black", linewidth=1)
    plt.tight_layout()
    plt.savefig(f"{feature_1}_{feature_2}_precursors_chromatogram.png", dpi=300)   
    plt.clf()  
    plt.close('all')
    plt.close('all')  


def mobilogram_overlay_precursors_nosmoothing(feature_1, feature_2, feature_df, tims_data):
    mobilogram_df = feature_df.loc[feature_df["feature_id"].isin([feature_1, feature_2])]
    im_array = tims_data.mobility_values

    plt.figure(figsize=(10,8))
    plt.style.use("seaborn")

    precursor_mz = mobilogram_df[mobilogram_df["feature_id"]==feature_1]["precursor_mz"].unique().item()
    
    RT_1 = mobilogram_df[mobilogram_df["feature_id"]==feature_1]["feature_retention_time"].unique().item()
    RT_2 = mobilogram_df[mobilogram_df["feature_id"]==feature_2]["feature_retention_time"].unique().item()

    IM_1 = mobilogram_df[mobilogram_df["feature_id"]==feature_1]["library_IM"].unique().item()
    IM_2 = mobilogram_df[mobilogram_df["feature_id"]==feature_2]["library_IM"].unique().item()

    peptide_1 = mobilogram_df[mobilogram_df["feature_id"]==feature_1]["FullPeptideName"].unique().item()
    peptide_2 = mobilogram_df[mobilogram_df["feature_id"]==feature_2]["FullPeptideName"].unique().item()
   
    feature_indices = get_feature_indices_overlay(tims_data, precursor_mz, IM_1, IM_2, RT_1, RT_2, im_tolerance=0.1)  
    feature_intensities = tims_data.bin_intensities(feature_indices, ["mobility_values"])

    non_zeros = np.flatnonzero(feature_intensities)
    start = max(0, non_zeros[0] - 1)
    end = non_zeros[-1] + 2
    im_values = im_array[start: end]
    feature_intensities = feature_intensities[start: end]

    plt.plot(im_values, feature_intensities, color="black", linestyle="-", alpha=1)

    feature_list = np.array([feature_1, feature_2])
    phosphorylation_positions = np.array([phosphorylation_position(peptide_1), phosphorylation_position(peptide_2)])
    phosphorylation_positions = np.argsort(phosphorylation_positions)
    features_sorted = feature_list[phosphorylation_positions]

    features_sorted = np.flip(features_sorted)

    for i in range(2):
        feature = features_sorted[i]
        IM = mobilogram_df[mobilogram_df["feature_id"]==feature]["library_IM"].unique().item()
        peptide = mobilogram_df[mobilogram_df["feature_id"]==feature]["FullPeptideName"].unique().item()
        
        if i ==1:
            rt_color = "deepskyblue" 
            plt.axvline(x=IM, color=rt_color, linewidth=1.5, linestyle=":", label=f"{peptide}")
            #plt.text(IM_1, np.amax(feature_intensities) + 100, f"Library IM = {IM}", horizontalalignment='left')
        else:     
            rt_color = "forestgreen"
            plt.axvline(x=IM, color=rt_color, linewidth=1.5, linestyle=":", label=f"{peptide}")
            #plt.text(IM, np.amax(feature_intensities) - np.amax(feature_intensities)*0.05, f"Library IM = {IM}", horizontalalignment='right')

    plt.legend(bbox_to_anchor=(1.05,1), loc="upper left", fontsize=10)
    plt.xlabel("Ion Mobility Vs $cm^{-2}$", fontsize=14)
    plt.ylabel("Intensity", fontsize=14)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    #plt.title(f"{peptide_1} and {peptide_2}", loc="left", fontsize=10)
    plt.axhline(y=0, color = "black", linewidth=1)
    plt.tight_layout()
    plt.savefig(f"{feature_1}_{feature_2}_raw_precursors_mobilograms.png", dpi=300)   
    plt.clf()  
    plt.close('all')


def mobilogram_precursor_phosphorylated_transitions_nosmoothing(feature, feature_df, tims_data):
 
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
        plt.plot(im_values, feature_intensities, color='black', linestyle="-", alpha=1, label=f"{peptide}")

    for i in range(len(mobilogram_df)): 
            transition = mobilogram_df.iloc[i]["transition_id"].item()
            precursor_mz = mobilogram_df[mobilogram_df["transition_id"]==transition]["precursor_mz"].item()
            transition_mz = mobilogram_df[mobilogram_df["transition_id"]==transition]["transition_mz"].item()
            rt = mobilogram_df[mobilogram_df["transition_id"]==transition]["feature_retention_time"].item()
            IM = mobilogram_df[mobilogram_df["transition_id"]==transition]["feature_IM"].item()
            transition_annotation = mobilogram_df[mobilogram_df["transition_id"]==transition]["transition_annotation"].item()
            
            transition_indices = get_transition_indices(tims_data, transition_mz, IM, rt, precursor_mz, im_tolerance=0.1)
            transition_intensities = tims_data.bin_intensities(transition_indices, ["mobility_values"])

            non_zeros = np.flatnonzero(transition_intensities)
            if len(non_zeros) <= 10:
                continue

            else:
                start = max(0, non_zeros[0] - 1)
                end = non_zeros[-1] + 2
                im_values = im_array[start: end]
                transition_intensities = transition_intensities[start: end]
    
                plt.plot(im_values, transition_intensities, color=clrs[i], linestyle="--", alpha=1, label=f"{transition_annotation}")
        

    plt.legend(bbox_to_anchor=(1.05,1), loc="upper left", fontsize=7)
    plt.xlabel("Ion Mobility Vs $cm^{-2}$")
    plt.ylabel("Intensity")
    plt.axhline(y=0, color = "black", linewidth=1)
    plt.tight_layout()
    plt.title(f"{peptide}", loc="left", fontsize=10)
    plt.savefig(f"{feature}_raw_mobilogram_precursor_phosphorylated_fragments.png", dpi=300, bbox_inches="tight")   
    plt.clf()
    plt.close('all')



def chromatogram_overlay_unique_phosphorylated_transitions_annotated(feature_1, feature_2, feature_df, tims_data, feature_1_pep, feature_2_pep, retention_time_delta):

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



    rt_array = tims_data.rt_values

    plt.figure(figsize=(10,8))
    plt.style.use("seaborn")

    plot_lines_1 = []
    plot_lines_2 = []
    feature_label_1 = []
    feature_label_2 = []

    plt.figure(figsize=(10,8))

    RT_1 = peptide_df.loc[peptide_df["feature_id"] == feature_1]["feature_retention_time"].unique().item()
    RT_2 = peptide_df.loc[peptide_df["feature_id"] == feature_2]["feature_retention_time"].unique().item()

    feature_list = np.array([feature_1, feature_2])
    phosphorylation_positions = np.array([phosphorylation_position(peptide_1), phosphorylation_position(peptide_2)])
    phosphorylation_positions = np.argsort(phosphorylation_positions)
    features_sorted = feature_list[phosphorylation_positions]

    features_sorted = np.flip(features_sorted)

    for j in range(2):
        feature = features_sorted[j]
        if j ==1: 
            clrs = sns.color_palette('Blues', n_colors= len(peptide_df.loc[peptide_df["feature_id"] == feature])+1)
        else:
            clrs = sns.color_palette('Greens', n_colors= len(peptide_df.loc[peptide_df["feature_id"] == feature])+1)
    
  
        precursor_mz = feature_df.loc[feature_df["feature_id"]==feature]["precursor_mz"].unique().item()
        transition_df_peptide = peptide_df.loc[peptide_df["feature_id"] == feature] 

        

        for i in range(len(transition_df_peptide)):
            if len(transition_df_peptide) > 0:
                transition = transition_df_peptide.iloc[i]["transition_id"]
                transition_mz = transition_df_peptide[transition_df_peptide["transition_id"]==transition]["transition_mz"].item()
                transition_annotation = transition_df_peptide[transition_df_peptide["transition_id"]==transition]["transition_annotation"].item()

                transition_indices = get_transition_indices_overlay(tims_data, transition_mz, 1,1, RT_1, RT_2, precursor_mz, im_tolerance=0.6, rt_tolerance=25)
                transition_intensities = tims_data.bin_intensities(transition_indices, ["rt_values"])
                

                non_zeros = np.flatnonzero(transition_intensities)
                if len(non_zeros) <= 10:
                    continue

                else:
                    start = max(0, non_zeros[0] - 1)
                    end = non_zeros[-1] + 2
                    rt_values = rt_array[start: end]
                    transition_intensities = transition_intensities[start: end]
                    smoothed_transition_intensities = savgol_filter(transition_intensities,11,3)

                    if feature == feature_1: 
                        l1, = plt.plot(rt_values, smoothed_transition_intensities, color=clrs[i+1], marker = ".", linestyle="-", alpha=1)
                        
                        plot_lines_1.append(l1)
                        feature_label_1.append(str(transition_annotation))
                    else: 
                        l2, = plt.plot(rt_values, smoothed_transition_intensities, color=clrs[i+1], marker = ".", linestyle="-", alpha=1)
                        plot_lines_2.append(l2)
                        feature_label_2.append(str(transition_annotation))
            else:
                break 

    plt.style.use("seaborn")
    legend_1 = plt.legend(plot_lines_1, feature_label_1, title=f"{peptide_1}", bbox_to_anchor=(1.05,1), loc="upper left",fontsize=7)
    legend_2 = plt.legend(plot_lines_2, feature_label_2, title=f"{peptide_2}", bbox_to_anchor=(1.05,0.35), loc="center left",fontsize=7)
    plt.gca().add_artist(legend_1)
    plt.gca().add_artist(legend_2)
    plt.xlabel("Retention time (s)")
    plt.ylabel("Intensity")
    plt.title(f"Retention time difference = {retention_time_delta}s \n IPF PEP scores: {peptide_1} = {feature_1_pep}, {peptide_2} = {feature_2_pep}", loc="left",fontsize=10)
    plt.axhline(y=0, color = "black", linewidth=1)
    plt.tight_layout()
    plt.savefig(f"{feature_1}_{feature_2}_smoothed_unique_phospho_fragments_chromatogram.png", dpi=300, bbox_inches="tight")   
    plt.clf()    
    plt.close('all')


def chromatogram_overlay_unique_phosphorylated_transitions_annotated_nosmoothing(feature_1, feature_2, feature_df, tims_data, feature_1_pep, feature_2_pep, retention_time_delta):


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
    rt_array = tims_data.rt_values

    plt.figure(figsize=(10,8))
    plt.style.use("seaborn")

    plot_lines_1 = []
    plot_lines_2 = []
    feature_label_1 = []
    feature_label_2 = []

    plt.figure(figsize=(10,8))

    RT_1 = peptide_df.loc[peptide_df["feature_id"] == feature_1]["feature_retention_time"].unique().item()
    RT_2 = peptide_df.loc[peptide_df["feature_id"] == feature_2]["feature_retention_time"].unique().item()

    feature_list = np.array([feature_1, feature_2])
    phosphorylation_positions = np.array([phosphorylation_position(peptide_1), phosphorylation_position(peptide_2)])
    phosphorylation_positions = np.argsort(phosphorylation_positions)
    features_sorted = feature_list[phosphorylation_positions]
    features_sorted = np.flip(features_sorted)

    for j in range(2):
        feature = features_sorted[j]
        if j ==1: 
            clrs = sns.color_palette('Blues', n_colors= len(peptide_df.loc[peptide_df["feature_id"] == feature])+1)
            rt_color = "deepskyblue"
        else:
            clrs = sns.color_palette('Greens', n_colors= len(peptide_df.loc[peptide_df["feature_id"] == feature])+1)
            rt_color = "forestgreen"
    
        precursor_mz = feature_df.loc[feature_df["feature_id"]==feature]["precursor_mz"].unique().item()
    
        transition_df_peptide = peptide_df.loc[peptide_df["feature_id"] == feature] 

        for i in range(len(transition_df_peptide)):
            if len(transition_df_peptide) > 0:
                transition = transition_df_peptide.iloc[i]["transition_id"]
                transition_mz = transition_df_peptide[transition_df_peptide["transition_id"]==transition]["transition_mz"].item()
                transition_annotation = transition_df_peptide[transition_df_peptide["transition_id"]==transition]["transition_annotation"].item()

                transition_indices = get_transition_indices_overlay(tims_data, transition_mz, 1,1, RT_1, RT_2, precursor_mz, im_tolerance=0.6)
                transition_intensities = tims_data.bin_intensities(transition_indices, ["rt_values"])
                

                non_zeros = np.flatnonzero(transition_intensities)
                if len(non_zeros) <= 10:
                    continue

                else:
                    start = max(0, non_zeros[0] - 1)
                    end = non_zeros[-1] + 2
                    rt_values = rt_array[start: end]
                    transition_intensities = transition_intensities[start: end]
                    transition_intensities[transition_intensities==0] = np.nan
                    transition_mask = np.isfinite(transition_intensities)

                    if feature == feature_1: 
                        l1, = plt.plot(rt_values[transition_mask], transition_intensities[transition_mask], color=clrs[i+1], marker = ".", linestyle="-", alpha=1)
                        plot_lines_1.append(l1)
                        feature_label_1.append(str(transition_annotation))
                        plt.axvline(x=RT_1, color=rt_color, linewidth=1, linestyle=":")
                    else: 
                        l2, = plt.plot(rt_values[transition_mask], transition_intensities[transition_mask], color=clrs[i+1], marker = ".", linestyle="-", alpha=1)
                        plot_lines_2.append(l2)
                        feature_label_2.append(str(transition_annotation))
                        plt.axvline(x=RT_2, color=rt_color, linewidth=1, linestyle=":")
            else:
                break 


    legend_1 = plt.legend(plot_lines_1, feature_label_1, title=f"{peptide_1}", bbox_to_anchor=(1.05,1), loc="upper left",fontsize=7)
    legend_2 = plt.legend(plot_lines_2, feature_label_2, title=f"{peptide_2}", bbox_to_anchor=(1.05,0.35), loc="center left",fontsize=7)
    plt.gca().add_artist(legend_1)
    plt.gca().add_artist(legend_2)
    plt.xlabel("Retention time (s)")
    plt.ylabel("Intensity")
    plt.title(f"Retention time difference = {retention_time_delta}s \n IPF PEP scores: {peptide_1} = {feature_1_pep}, {peptide_2} = {feature_2_pep}", loc="left",fontsize=10)
    plt.axhline(y=0, color = "black", linewidth=1)
    plt.tight_layout()
    plt.savefig(f"{feature_1}_{feature_2}_raw_unique_phospho_fragments_chromatogram.png", dpi=300, bbox_inches="tight")   
    plt.clf()    
    plt.close('all')


def chromatogram_overlay_precursors_unique_phosphorylated_transitions_annotated_nosmoothing(feature_1, feature_2, feature_df, tims_data, feature_1_pep, feature_2_pep, retention_time_delta):

    peptide_df_1 = feature_df.loc[feature_df["feature_id"] == feature_1]
    peptide_1 = feature_df.loc[feature_df["feature_id"] == feature_1]["FullPeptideName"].unique().item()
    peptide_df_1 = peptide_df_1.loc[peptide_df_1["transition_annotation_nocharge"].isin(fragments_with_phosphorylation(peptide_1))]

    peptide_df_2 = feature_df.loc[feature_df["feature_id"] == feature_2]
    peptide_2 = feature_df.loc[feature_df["feature_id"] == feature_2]["FullPeptideName"].unique().item()
    peptide_df_2 = peptide_df_2.loc[peptide_df_2["transition_annotation_nocharge"].isin(fragments_with_phosphorylation(peptide_2))]

    peptide_df = pd.concat([peptide_df_1, peptide_df_2], axis=0)
    peptide_df = peptide_df.loc[peptide_df["transition_identifying"]==1]

    peptide_df["n_features"] = (peptide_df.groupby(["transition_annotation_nocharge"])["feature_id"].transform("nunique")) 
    peptide_df = peptide_df.loc[peptide_df["n_features"] != 2] 
    #peptide_df = peptide_df.loc[peptide_df.groupby(["feature_id", "transition_annotation", "transition_mz"])["transition_identifying"].idxmax()]

    rt_array = tims_data.rt_values

    plt.figure(figsize=(10,8))
    plt.style.use("seaborn")

    RT_1 = peptide_df.loc[peptide_df["feature_id"] == feature_1]["feature_retention_time"].unique().item()
    RT_2 = peptide_df.loc[peptide_df["feature_id"] == feature_2]["feature_retention_time"].unique().item()
    IM_1 = peptide_df.loc[peptide_df["feature_id"] == feature_1]["feature_IM"].unique().item()
    IM_2 = peptide_df.loc[peptide_df["feature_id"] == feature_2]["feature_IM"].unique().item()

    precursor_mz = feature_df.loc[feature_df["feature_id"]==feature_1]["precursor_mz"].unique().item()
 
    feature_indices = get_feature_indices_overlay(tims_data, precursor_mz, 1, 1, RT_1, RT_2, rt_tolerance=60, im_tolerance=0.6)       # don't want to consider ion mobility while doing chromatogram extraction
    
    x = tims_data.as_dataframe(feature_indices)
    x = x.loc[(x["frame_indices"]-1) % 17 ==0]
    feature_indices = x["raw_indices"].unique()

    feature_intensities = tims_data.bin_intensities(feature_indices, ["rt_values"])

    non_zeros = np.flatnonzero(feature_intensities)
    start = max(0, non_zeros[0] - 1)
    end = non_zeros[-1] + 2
    rt_values = rt_array[start: end]
    feature_intensities = feature_intensities[start: end]

    feature_intensities[feature_intensities==0] = np.nan
    feature_mask = np.isfinite(feature_intensities)

    plot_lines_1 = []
    plot_lines_2 = []
    feature_label_1 = []
    feature_label_2 = []

    l0 = plt.plot(rt_values[feature_mask], feature_intensities[feature_mask], color="black", linestyle="-", alpha=1)      
    plot_lines_1.append(l0)
    feature_label_1.append("precursor")
    plot_lines_2.append(l0)
    feature_label_2.append("precursor")

    feature_list = np.array([feature_1, feature_2])
    phosphorylation_positions = np.array([phosphorylation_position(peptide_1), phosphorylation_position(peptide_2)])
    phosphorylation_positions = np.argsort(phosphorylation_positions)
    features_sorted = feature_list[phosphorylation_positions]
    features_sorted = np.flip(features_sorted)


    for j in range(2):
        feature = features_sorted[j]
        peptide = feature_df.loc[feature_df["feature_id"]==feature]["FullPeptideName"].unique().item()
        precursor_mz = feature_df.loc[feature_df["feature_id"]==feature]["precursor_mz"].unique().item()
        transition_df_peptide = peptide_df.loc[peptide_df["feature_id"] == feature] 
        RT = peptide_df.loc[peptide_df["feature_id"] == feature]["feature_retention_time"].unique().item()
        IM = peptide_df.loc[peptide_df["feature_id"] == feature]["feature_IM"].unique().item()
        if j ==1: 
            clrs = sns.color_palette('Blues', n_colors= len(transition_df_peptide)+1)
            rt_color = "deepskyblue"
            plt.axvline(x=RT, color=rt_color, linewidth=1.5, linestyle=":")

            for i in range(len(transition_df_peptide)):
                if len(transition_df_peptide) > 0:
                    transition = transition_df_peptide.iloc[i]["transition_id"]
                    transition_mz = transition_df_peptide[transition_df_peptide["transition_id"]==transition]["transition_mz"].item()
                    transition_annotation = transition_df_peptide[transition_df_peptide["transition_id"]==transition]["transition_annotation"].item()

                    transition_indices = get_transition_indices_overlay(tims_data, transition_mz, IM_1,IM_2, RT_1, RT_2, precursor_mz, rt_tolerance=60, im_tolerance=0.05)
                    transition_intensities = tims_data.bin_intensities(transition_indices, ["rt_values"])

                    non_zeros = np.flatnonzero(transition_intensities)
                    if len(non_zeros) <= 10:
                        continue

                    else:
                        start = max(0, non_zeros[0] - 1)
                        end = non_zeros[-1] + 2
                        rt_values = rt_array[start: end]
                        transition_intensities = transition_intensities[start: end]
                        transition_intensities[transition_intensities==0] = np.nan
                        transition_mask = np.isfinite(transition_intensities)

                        l2, = plt.plot(rt_values[transition_mask], transition_intensities[transition_mask], color=clrs[i+1], linestyle="-", alpha=1)
                        plot_lines_2.append(l2)
                        feature_label_2.append(str(transition_annotation))
                else:
                    break 
            
            legend_2 = plt.legend(plot_lines_2, feature_label_2, title=f"{peptide}", bbox_to_anchor=(1.05,0.20), loc="center left",fontsize=7)

        
        else:
            clrs = sns.color_palette('Greens', n_colors= len(transition_df_peptide)+1)
            rt_color = "forestgreen"
            plt.axvline(x=RT, color=rt_color, linewidth=1.5, linestyle=":")

            for i in range(len(transition_df_peptide)):
                if len(transition_df_peptide) > 0:
                    transition = transition_df_peptide.iloc[i]["transition_id"]
                    transition_mz = transition_df_peptide[transition_df_peptide["transition_id"]==transition]["transition_mz"].item()
                    transition_annotation = transition_df_peptide[transition_df_peptide["transition_id"]==transition]["transition_annotation"].item()

                    transition_indices = get_transition_indices_overlay(tims_data, transition_mz, IM, IM, RT_1, RT_2, precursor_mz, rt_tolerance=60, im_tolerance=0.6)
                    transition_intensities = tims_data.bin_intensities(transition_indices, ["rt_values"])

                    non_zeros = np.flatnonzero(transition_intensities)
                    if len(non_zeros) <= 10:
                        continue

                    else:
                        start = max(0, non_zeros[0] - 1)
                        end = non_zeros[-1] + 2
                        rt_values = rt_array[start: end]
                        transition_intensities = transition_intensities[start: end]
                        transition_intensities[transition_intensities==0] = np.nan
                        transition_mask = np.isfinite(transition_intensities)

                        l1, = plt.plot(rt_values[transition_mask], transition_intensities[transition_mask], color=clrs[i+1], linestyle="-", alpha=1)
                        plot_lines_1.append(l1)
                        feature_label_1.append(str(transition_annotation))
                else:
                    break 
            
            legend_1 = plt.legend(plot_lines_1, feature_label_1, title=f"{peptide}", bbox_to_anchor=(1.05,1), loc="upper left",fontsize=7)

        
    plt.gca().add_artist(legend_1)
    plt.gca().add_artist(legend_2)
    plt.xlabel("Retention time (s)", fontsize=14)
    plt.ylabel("Intensity", fontsize=14)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    #plt.title(f"{peptide_1} and {peptide_2}", loc="left",fontsize=10)
    plt.axhline(y=0, color = "black", linewidth=1)
    plt.tight_layout()
    plt.savefig(f"{feature_1}_{feature_2}_raw_unique_phospho_fragments_precursors_chromatogram.png", dpi=300, bbox_inches="tight")   
    plt.clf()    
    plt.close('all')

def mobilogram_noncoeluting_overlay_precursors_unique_phosphorylated_transitions_annotated_nosmoothing(feature_1, feature_2, feature_df, tims_data):

    peptide_df_1 = feature_df.loc[feature_df["feature_id"] == feature_1]
    peptide_1 = feature_df.loc[feature_df["feature_id"] == feature_1]["FullPeptideName"].unique().item()
    peptide_df_1 = peptide_df_1.loc[peptide_df_1["transition_annotation_nocharge"].isin(fragments_with_phosphorylation(peptide_1))]

    peptide_df_2 = feature_df.loc[feature_df["feature_id"] == feature_2]
    peptide_2 = feature_df.loc[feature_df["feature_id"] == feature_2]["FullPeptideName"].unique().item()
    peptide_df_2 = peptide_df_2.loc[peptide_df_2["transition_annotation_nocharge"].isin(fragments_with_phosphorylation(peptide_2))]

    peptide_df = pd.concat([peptide_df_1, peptide_df_2], axis=0)
    peptide_df = peptide_df.loc[peptide_df["transition_identifying"]==1]
    peptide_df["n_features"] = (peptide_df.groupby(["transition_annotation_nocharge"])["feature_id"].transform("nunique")) 
    peptide_df = peptide_df.loc[peptide_df["n_features"] != 2] 
    #peptide_df = peptide_df.loc[peptide_df.groupby(["feature_id", "transition_annotation", "transition_mz"])["transition_identifying"].idxmax()]

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
    features_sorted = np.flip(features_sorted)

    for j in range(2):
        feature = features_sorted[j]
        if j ==1: 
            clrs = sns.color_palette('Blues', n_colors= len(peptide_df.loc[peptide_df["feature_id"] == feature])+1)
        else:
            clrs = sns.color_palette('Greens', n_colors= len(peptide_df.loc[peptide_df["feature_id"] == feature])+1)
    
        precursor_mz = feature_df.loc[feature_df["feature_id"]==feature]["precursor_mz"].unique().item()
        RT = peptide_df.loc[peptide_df["feature_id"] == feature]["feature_retention_time"].unique().item()
        IM = peptide_df.loc[peptide_df["feature_id"] == feature]["library_IM"].unique().item()

        transition_df_peptide = peptide_df.loc[peptide_df["feature_id"] == feature] 

        for i in range(len(transition_df_peptide)):
            if len(transition_df_peptide) > 0:
                transition = transition_df_peptide.iloc[i]["transition_id"]
                transition_mz = transition_df_peptide[transition_df_peptide["transition_id"]==transition]["transition_mz"].item()
                transition_annotation = transition_df_peptide[transition_df_peptide["transition_id"]==transition]["transition_annotation"].item()

                transition_indices = get_transition_indices(tims_data, transition_mz, IM, RT, precursor_mz, im_tolerance=0.15)
                transition_intensities = tims_data.bin_intensities(transition_indices, ["mobility_values"])

                non_zeros = np.flatnonzero(transition_intensities)
                if len(non_zeros) <= 10:
                    continue

                else:
                    start = max(0, non_zeros[0] - 1)
                    end = non_zeros[-1] + 2
                    im_values = im_array[start: end]
                    transition_intensities = transition_intensities[start: end]

                    if feature == feature_1: 
                        l1, = plt.plot(im_values, transition_intensities, color=clrs[i+1], linestyle="--", alpha=1)
                        plot_lines_1.append(l1)
                        feature_label_1.append(str(transition_annotation))
                    else: 
                        l2, = plt.plot(im_values, transition_intensities, color=clrs[i+1], linestyle="--", alpha=1)
                        plot_lines_2.append(l2)
                        feature_label_2.append(str(transition_annotation))
            else:
                break 
        
        feature_indices = get_feature_indices(tims_data,precursor_mz,  IM, RT, im_tolerance=0.15)
        feature_intensities = tims_data.bin_intensities(feature_indices, ["mobility_values"])

        non_zeros = np.flatnonzero(feature_intensities)
        if len(non_zeros) == 0:
            continue
        else:
            start = max(0, non_zeros[0] - 1)
            end = non_zeros[-1] + 2
            im_values = im_array[start: end]
            feature_intensities = feature_intensities[start: end]
            feature_intensities = feature_intensities*0.1

            if feature == feature_1: 
                l1, = plt.plot(im_values, feature_intensities, color="black", linestyle="-", alpha=0.8)
                #plt.text(IM, np.amax(feature_intensities) + 100, f"Library IM = {IM}", horizontalalignment='center')
                plot_lines_1.append(l1)
                feature_label_1.append("precursor")

            else: 
                l2, = plt.plot(im_values, feature_intensities, color="black", linestyle="-", alpha=0.8)
                #plt.text(IM, np.amax(feature_intensities) + 100, f"Library IM = {IM}", horizontalalignment='center')
                plot_lines_2.append(l2)
                feature_label_2.append("precursor")



    plt.style.use("seaborn")
    legend_1 = plt.legend(plot_lines_1, feature_label_1, title=f"{peptide_1}", bbox_to_anchor=(1.05,1), loc="upper left",fontsize=7)
    legend_2 = plt.legend(plot_lines_2, feature_label_2, title=f"{peptide_2}", bbox_to_anchor=(1.05,0.20), loc="center left",fontsize=7)
    plt.gca().add_artist(legend_1)
    plt.gca().add_artist(legend_2)
    plt.xlabel("Ion Mobility Vs $cm^{-2}$", fontsize=14)
    plt.ylabel("Intensity", fontsize=14)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    #plt.title(f"Retention time difference = {retention_time_delta}s \n IPF PEP scores: {peptide_1} = {feature_1_pep}, {peptide_2} = {feature_2_pep}", loc="left",fontsize=10)
    plt.axhline(y=0, color = "black", linewidth=1)
    plt.tight_layout()
    plt.savefig(f"{feature_1}_{feature_2}_noncoeluting_raw_unique_phospho_fragments_precursors_mobilogram.png", dpi=300, bbox_inches="tight")   
    plt.clf()    
    plt.close('all')