#!/bin/bash

library=upool_library.pqp
irt=hroest_DIA_iRT.TraML

for file in /scratch/h/hroest/asrini/timstofdata_mzml/*.mzML; 
do mzml=$(basename $file);
sbatch run_openswathworkflow_1.sh $mzml $library $irt; 
done