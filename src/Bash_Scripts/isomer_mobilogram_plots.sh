ISOMER_DATAFRAME=
PYPROPHET_FOLDER=

source /home/roestlab/anaconda3/etc/profile.d/conda.sh 
conda activate alphatims

for i in /media/roestlab/Data1/User/AparnaS/data/TimsTOF/M_pool_diaPASEF/*d; 
python /media/roestlab/Data1/User/AparnaS/Srinivasan_MS_Thesis/src/Synthetic_Phosphopeptide_Isomers_Analysis_Python/isomer_transition_data.py --path-to-raw-data $i --isomer-df $ISOMER_DATAFRAME --data-folder $PYPROPHET_FOLDER
done; 

conda deactivate 