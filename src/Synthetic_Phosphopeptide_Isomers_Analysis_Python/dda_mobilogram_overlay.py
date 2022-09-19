from alphatims_functions import *
from export_format_functions import * 
import click 

@click.command()
@click.option('--feature-1', help="Feature number one", type=int)
@click.option('--feature-2', help="Feature number two", type=int)
@click.option('--path-to-raw-data', help="Path to raw data file",  type=str, required=True)
@click.option('--feature-dataframe', help="Path to feature isomer dataframe",  type=click.Path(exists=True), required=True)
def main(feature_1, feature_2, feature_dataframe, path_to_raw_data):
    feature_dataframe = pd.read_csv(feature_dataframe, sep="\t")
    feature_df = feature_dataframe.loc[(feature_dataframe["feature_id.x"] ==feature_1) & (feature_dataframe["feature_id.y"]==feature_2)]

    filename_1 = feature_df["filename.x"].item()
    filename_2 = feature_df["filename.y"].item()

    tims_data_1 = alphatims.bruker.TimsTOF(f"{path_to_raw_data}/{filename_1}.d")
    tims_data_2 = alphatims.bruker.TimsTOF(f"{path_to_raw_data}/{filename_2}.d")

    plt.figure(figsize=(10,8))
    
    mz = feature_df["mz"].item()
    peptide_1 = feature_df["FullPeptideName.x"].item()
    RT_1 = feature_df["RT.x"].item()
    IM_1 = feature_df["IM.x"].item()
    library_IM_1 = feature_df["library_IM.x"].item()
    feature_indices_1 = get_feature_indices(tims_data_1, mz, IM_1, RT_1, rt_tolerance=20, im_tolerance = 0.1)

    peptide_2 = feature_df["FullPeptideName.y"].item()
    RT_2 = feature_df["RT.y"].item()
    IM_2 = feature_df["IM.y"].item()
    library_IM_2 = feature_df["library_IM.y"].item()
    feature_indices_2 = get_feature_indices(tims_data_2, mz, IM_2, RT_2, rt_tolerance=20, im_tolerance = 0.1)

    im_array = tims_data_1.mobility_values

    feature_list = np.array([feature_1, feature_2])
    phosphorylation_positions = np.array([phosphorylation_position(peptide_1), phosphorylation_position(peptide_2)])
    phosphorylation_positions = np.argsort(phosphorylation_positions)
    features_sorted = feature_list[phosphorylation_positions]
    clrs = ["deepskyblue", "forestgreen"]


    for i in range(2):
        feature = features_sorted[i]
        clr = clrs[i]
       
        if feature == feature_1:
            feature_intensities = tims_data_1.bin_intensities(feature_indices_1, ["mobility_values"]) 
            non_zeros = np.flatnonzero(feature_intensities)
            start = max(0, non_zeros[0] - 1)
            end = non_zeros[-1] + 2
            im_values = im_array[start: end]
            feature_intensities = feature_intensities[start: end]
            
            plt.plot(im_values, feature_intensities, color=clr,linestyle="-", alpha=1)
            plt.axvline(x=library_IM_1, color=clr, linewidth=1, linestyle=":", label=f"{peptide_1}")
            plt.text(library_IM_1, np.amax(feature_intensities) + np.amax(feature_intensities)*0.01, f"Library IM = {library_IM_1}", horizontalalignment='left')
    
        else: 
            feature_intensities = tims_data_2.bin_intensities(feature_indices_2, ["mobility_values"])
            non_zeros = np.flatnonzero(feature_intensities)
            start = max(0, non_zeros[0] - 1)
            end = non_zeros[-1] + 2
            im_values = im_array[start: end]
            feature_intensities = feature_intensities[start: end]

            plt.plot(im_values, feature_intensities, color=clr, linestyle="-", alpha=1)
            plt.axvline(x=library_IM_2, color=clr, linewidth=1, linestyle=":", label=f"{peptide_2}")
            plt.text(library_IM_2, np.amax(feature_intensities) - np.amax(feature_intensities)*0.01, f"Library IM = {library_IM_2}", horizontalalignment='right')
        
    
    plt.style.use("seaborn")
    plt.legend(bbox_to_anchor=(1.05,1), loc="upper left", fontsize=7)
    plt.xlabel("Ion Mobility Vs $cm^{-2}$")
    plt.ylabel("Intensity")
    plt.title(f"{peptide_1} and {peptide_2}", loc="left", fontsize=9)
    plt.axhline(y=0, color = "black", linewidth=1)
    plt.tight_layout()
    plt.savefig(f"{feature_1}_{feature_2}_dda_raw_precursors_mobilograms.png", dpi=300)   
    plt.clf()  
    plt.close('all')



if __name__ == "__main__":
    main()