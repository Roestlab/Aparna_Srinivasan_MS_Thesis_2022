source /home/roestlab/anaconda3/etc/profile.d/conda.sh && 
conda activate aparna
for i in /home/roestlab/aparna_diapasef/results/openswath/*.osw;
do file=$(basename $i);
tmp_dir=$(cat /proc/sys/kernel/random/uuid) && mkdir $tmp_dir; 
cp /home/roestlab/aparna_diapasef/results/openswath/$file $tmp_dir/; 
pyprophet score --in $tmp_dir/$file --out $tmp_dir/$file --classifier=LDA --level=ms2 --xeval_num_iter=3 --ss_initial_fdr=0.2 --ss_iteration_fdr=0.01 --threads=4;
pyprophet score --in $tmp_dir/$file --classifier=LDA --level=transition --ipf_min_transition_sn=-1 --xeval_num_iter=3 --ss_initial_fdr=0.2 --ss_iteration_fdr=0.01 --threads=4;
pyprophet ipf --no-ipf_ms1_scoring --no-ipf_ms2_scoring --ipf_grouped_fdr --in $tmp_dir/$file;
mv $tmp_dir/* /home/roestlab/aparna_diapasef/results/;
#mv $tmp_dir/*.osw /home/roestlab/aparna_diapasef/results/pyprophet/runspecific;
#mv $tmp_dir/*.pdf /home/roestlab/aparna_diapasef/results/pyprophet/runspecific;
rm -r $tmp_dir;
done
