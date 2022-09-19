source activate /home/roestlab/anaconda3/envs/josh_maxquant/
mono --aot=full /home/roestlab/anaconda3/envs/josh_maxquant/bin/MaxQuantCmd.exe
mono /home/roestlab/anaconda3/envs/josh_maxquant/bin/MaxQuantCmd.exe mqpar_upools_concatenatedfasta.xml
conda deactivate
 
