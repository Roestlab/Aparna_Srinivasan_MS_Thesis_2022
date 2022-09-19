mpool_list=("M1" "M2" "M3" "M4" "M5" "M6" "M7" "M8" "M9" "M10" "M11" "M12" "M13" "M14" "M15" "M16")
unimodfile=/home/roestlab/aparna_diapasef/unimod_ipf.xml
path=/home/roestlab/Development/openms/build/openms_v25/bin
data_dir=/media/roestlab/Data1/User/AparnaS/Synthetic_phosphopeptide_isomers_TimsTOF/Upools_MQ_speclib/mpool_specific_libraries
for i in ${mpool_list[@]};
do transition_list=${i}_transition_list.tsv; 
cp $data_dir/$transition_list /home/roestlab/aparna_diapasef/;
$path/OpenSwathAssayGenerator -in /home/roestlab/aparna_diapasef/$transition_list -out /home/roestlab/aparna_diapasef/${i}_lib_targets.pqp -min_transitions 6 -max_transitions 6 -product_lower_mz_limit 300 -product_upper_mz_limit 1800 -enable_ipf -unimod_file $unimodfile -disable_identification_ms2_precursors -disable_identification_specific_losses && $path/OpenSwathDecoyGenerator -in /home/roestlab/aparna_diapasef/${i}_lib_targets.pqp -out /home/roestlab/aparna_diapasef/${i}_library.pqp;
done 
