#!/bin/bash

# Author: Alex Weber, Ph.D.
# Date: Dec 11 2020


bold=mean_func
t1=../../../T1w/T1w_acpc_dc_restore_1.25


bet $t1 ${t1}_brain -R -m

epi_reg --epi=${bold}_brain.nii.gz --t1=${t1}.nii.gz --t1brain=${t1}_brain.nii.gz --out=BOLD-to-T1 #registration BOLD to T1; Time: 3mins

flirt -ref ${FSLDIR}/data/standard/MNI152_T1_2mm_brain -in ${t1}_brain.nii.gz -out T1-to-MNI -omat T1-to-MNI.mat #registration T1 to ref; Time:

convert_xfm -omat BOLD-to-MNI.mat -concat T1-to-MNI.mat BOLD-to-T1.mat

flirt -in filtered_func_data_clea_Hurst_Welch.nii.gz -ref ${FSLDIR}/data/standard/MNI152_T1_2mm_brain -out Welch-to-MNI -init BOLD-to-MNI.mat -applyxfm #register BOLD to MNI space

slicer Welch-to-MNI.nii.gz /usr/local/fsl/data/standard/MNI152_T1_2mm_brain.nii.gz -a reg.ppm
