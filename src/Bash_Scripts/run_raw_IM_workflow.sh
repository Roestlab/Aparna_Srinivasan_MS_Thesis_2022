#!/bin/bash
source /home/roestlab/anaconda3/etc/profile.d/conda.sh

conda activate alphatims

dataframe=/media/roestlab/Data1/User/AparnaS/Srinivasan_MS_Thesis/Data/Data_Analysis_Results/filtered_diapasef_results_5_1.tsv

for i in /media/roestlab/Data1/User/AparnaS/data/TimsTOF/M_pool_diaPASEF/*.d; 
do filename=$(basename $i .d); 
output_folder = /media/roestlab/Data1/User/AparnaS/Srinivasan_MS_Thesis/Data/Data_Analysis_Results/$filename; 
mkdir $output_folder;
python /media/roestlab/Data1/User/AparnaS/Srinivasan_MS_Thesis/src/Synthetic_Phosphopeptide_Isomers_Analysis_Python/raw_IM_extraction.py --path-to-raw-data $i --df-folder $dataframe --output-folder $output_folder; 
done 

conda deactivate