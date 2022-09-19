source /home/roestlab/anaconda3/etc/profile.d/conda.sh 
conda activate aparna 

path=/home/roestlab/anaconda3/envs/aparna/bin 

for i in /media/roestlab/Data1/User/AparnaS/data/TimsTOF/diapasef_100/*.d;
do tmp_dir=$(cat /proc/sys/kernel/random/uuid);
mkdir $tmp_dir; 
cp -r $i $tmp_dir/;
mzml_name=$(basename $i .d); 
$path/convertTDFtoMzML.py -a=$tmp_dir/$(basename $i) --output_name=$mzml_name.mzml --overlap=2;
mv *.mzML /media/roestlab/Data1/User/AparnaS/data/TimsTOF/diapasef_mzml_converted; 
rm -r $tmp_dir;
done

conda deactivate
