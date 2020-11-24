# FractalDimension

## BOLD_masks.sh 

## BOLD_masks_7T3T.sh

## FracTool_Call_3D.py 

FracTool_Call_3D.py runs the FracTool_Final.py script on every voxel, outputting the Hurst exponent for each voxel.

It is run with './FracTool_Call_3D.py filename.nii.gz'and the output is a filename_Hurst.nii.gz file, which is a 3D image of the Hurst exponent.  


## dfa_script_final.py

dfa_script_final.py runs the nolds function, DFA on every voxel in the brain and ignores empty voxels.

It is run with './dfa_script_final.py filename.nii.gz' and the output is a filename filename_Hurst_DFA.nii.gz file, which is a 3D image of the Hurst exponent.

## fmriprep_masks.sh

## hurst_rs_script_final.py

hurst_rs_script_final.py runs the nolds function, HurstR/S on every voxel in the brain.

It is run with './hurst_rs_script_final.py filename.nii.gz' and the output is a filename filename_Hurst_HurstRS.nii.gz file, which is a 3D image of the Hurst exponent.
