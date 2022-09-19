#!/bin/bash
#SBATCH --nodes=1
#SBATCH --account=def-hroest
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=4

module load cmake
module load autotools
module load gcc/8.3.0
module load qt/5.12.4

data_dir=/scratch/h/hroest/asrini/timstofdata_mzml/
analysis_dir=/home/h/hroest/asrini/data/Synthetic_phosphopeptide_isomers/timstof_mpools_analysis
path=/home/h/hroest/asrini/bin/openms/OpenMS-build-2.6.0/bin
 
mzml=$1
library=$2
Irtfile=$3

$path/OpenSwathWorkflow -in $SLURM_TMPDIR/$mzml -tr $SLURM_TMPDIR/$library -tr_irt $SLURM_TMPDIR/$Irtfile -enable_ms1 true -enable_ipf true -out_osw $filename.osw -out_chrom $filename.sqMass -rt_extraction_window 250 -mz_extraction_window 25 -mz_extraction_window_unit ppm -mz_extraction_window_ms1 25 -mz_extraction_window_ms1_unit ppm -use_ms1_ion_mobility false -irt_mz_extraction_window_unit ppm -irt_mz_extraction_window 40 -min_coverage 0.1 -readOptions cacheWorkingInMemory -mz_correction_function quadratic_regression_delta_ppm -tempDirectory $SLURM_TMPDIR -Debugging:irt_trafo $filename.trafoXML -log $filename.log -debug 10 -threads 4 -RTNormalization:alignmentMethod lowess -RTNormalization:NrRTBins 8 -RTNormalization:MinBinsFilled 4 -RTNormalization:lowess:span 0.01 -Scoring:stop_report_after_feature 5 -Scoring:TransitionGroupPicker:PeakPickerMRM:sgolay_frame_length 11 -Scoring:Scores:use_uis_scores
 
