
from alphatims_functions import *

def raw_ion_mobility_data_extraction(tims_data, im_extraction_df ,output_folder):

    for feature in im_extraction_df["feature_id"]:
        mz = float(im_extraction_df.loc[im_extraction_df["feature_id"] == feature]["mz"].item())
        IM = float(im_extraction_df.loc[im_extraction_df["feature_id"] == feature]["library_IM"].item())
        RT = float(im_extraction_df.loc[im_extraction_df["feature_id"] == feature]["RT"].item())    
         
        fullpeptidename =  im_extraction_df.loc[im_extraction_df["feature_id"] == feature]["FullPeptideName"].item()
        sequence = im_extraction_df.loc[im_extraction_df["feature_id"] == feature]["Sequence"].item()
        charge = str(im_extraction_df.loc[im_extraction_df["feature_id"] == feature]["Charge"].item())
        feature_indices = get_feature_indices(tims_data, mz, IM, RT)

        feature_im_values = tims_data.mobility_values 
        del_im = round(abs(np.mean(np.diff(feature_im_values))),6)
        feature_intensities = tims_data.bin_intensities(feature_indices, ["mobility_values"])   # this is a numpy array this bins across ALL tims scan values , so the size will always be the same as mobility_values
        non_zeros = np.flatnonzero(feature_intensities)

        if len(non_zeros) == 0:  # no non 0 intensity values 
            feature_im_values = np.empty(0, dtype=x_ticks.dtype)  # modify 
            feature_intensities = np.empty(0, dtype=intensities.dtype) #modify 

            im_extraction_df.loc[im_extraction_df["feature_id"] == feature, "n_peaks"] = np.nan
            im_extraction_df.loc[im_extraction_df["feature_id"] == feature, "smoothed_IM_peak"] = np.nan
            im_extraction_df.loc[im_extraction_df["feature_id"] == feature, "smoothed_IM_peakwidth"] = np.nan
            im_extraction_df.loc[im_extraction_df["feature_id"] == feature, "smoothed_peak_intensity"] = np.nan
            im_extraction_df.loc[im_extraction_df["feature_id"] == feature, "raw_max_intensity"] = np.nan
            im_extraction_df.loc[im_extraction_df["feature_id"] == feature, "peak_prominence"] = np.nan

        else:
            start = max(0, non_zeros[0] - 1)
            end = non_zeros[-1] + 2
            feature_im_values = feature_im_values[start: end]
            feature_intensities = feature_intensities[start: end]
            im_extraction_df.loc[im_extraction_df["feature_id"] == feature, "raw_max_intensity"] = np.amax(feature_intensities)
            # SMOOTHING 
            smoothed_intensity = savgol_filter(feature_intensities, 25,2)

           
            if np.amax(smoothed_intensity) < 5000:
                peaks = find_peaks(smoothed_intensity, prominence=500)[0]
                im_extraction_df.loc[im_extraction_df["feature_id"] == feature, "n_peaks"] = len(peaks)
    
            else:
                peaks = find_peaks(smoothed_intensity, prominence=2000)[0]   
                im_extraction_df.loc[im_extraction_df["feature_id"] == feature, "n_peaks"] = len(peaks)
            
            if len(peaks) == 0:
                im_extraction_df.loc[im_extraction_df["feature_id"] == feature, "smoothed_IM_peak"] = np.nan
                im_extraction_df.loc[im_extraction_df["feature_id"] == feature, "smoothed_IM_peakwidth"] = np.nan
                im_extraction_df.loc[im_extraction_df["feature_id"] == feature, "smoothed_peak_intensity"] = np.nan
                im_extraction_df.loc[im_extraction_df["feature_id"] == feature, "peak_prominence"] = np.nan
            else:
                idx = np.argmin(abs(feature_im_values[peaks] - IM)) # when there are multiple peaks, the peak closest to library value is taken
                im_extraction_df.loc[im_extraction_df["feature_id"] == feature, "smoothed_peak_intensity"] = smoothed_intensity[peaks[idx]] 
                smoothed_im_peak = feature_im_values[peaks[idx]]
                im_extraction_df.loc[im_extraction_df["feature_id"] == feature, "smoothed_IM_peak"] = smoothed_im_peak 
                im_extraction_df.loc[im_extraction_df["feature_id"] == feature, "smoothed_IM_peakwidth"]  = float(peak_widths(smoothed_intensity, np.array([peaks[idx]]), rel_height=0.5)[0]*del_im)
                im_extraction_df.loc[im_extraction_df["feature_id"] == feature, "peak_prominence"] = float(peak_prominences(smoothed_intensity,  np.array([peaks[idx]]))[0])

                #plt.style.use("seaborn")  

                #plt.plot(feature_im_values, feature_intensities, label="Raw values")
                #plt.plot(feature_im_values, smoothed_intensity, label="Smoothed Values") 
                #plt.text(IM, smoothed_intensity[peaks[idx]] + 100, f"Library IM = {IM}", horizontalalignment='center')
                #plt.xlabel("Ion Mobility Vs $cm^{-2}$")
                #plt.ylabel("Intensity")
                #plt.title(f"{fullpeptidename}\nCharge state={charge}\nIM={smoothed_im_peak}", loc="left")
                #plt.figlegend()
                #plt.savefig(f"{output_folder}/{feature}_{sequence}.png", dpi=300)
                #plt.clf()
                #plt.close('all')


    return im_extraction_df

