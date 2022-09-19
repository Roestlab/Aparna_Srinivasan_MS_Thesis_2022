python /home/roestlab/aparna_diaPASEF/Synthetic_phosphopeptide_isomers_TimsTOF/evidence_file_processing.py \
--evidence-file /home/roestlab/aparna_diaPASEF/aparna_ddapasef/ddapasef_upool/combined/txt/evidence.txt \
--output-file /media/roestlab/Data1/User/AparnaS/Synthetic_phosphopeptide_isomers_TimsTOF/Upools_MQ_speclib/upool_library/evidence.txt \
--min-probability 0.75 \
--ground-truth /media/roestlab/Data1/User/AparnaS/Synthetic_phosphopeptide_isomers_TimsTOF/Upools_MQ_speclib/upool_library/ground_truth_table.tsv 

python /home/roestlab/aparna_diaPASEF/Synthetic_phosphopeptide_isomers_TimsTOF/reformat_msms.py \
--msms-file /home/roestlab/aparna_diaPASEF/aparna_ddapasef/ddapasef_upool/combined/txt/msms.txt \
--msms-output /media/roestlab/Data1/User/AparnaS/Synthetic_phosphopeptide_isomers_TimsTOF/Upools_MQ_speclib/upool_library/msms.txt \

