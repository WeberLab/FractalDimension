#!/bin/bash

# Author: Alex Weber, Ph.D.
# Date: Oct 9 2020

if [ $# -lt 1 ]
then
    echo "Usage: $0 <fmriprepBOLD brain>"
    echo "Example: $0 sub-01_task-rest_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz"
    echo "This grabs the T1 from MNI space, then creates very conservative whitematter and greymatter masks"
    echo "[note: default cluster removal is 10 and 6]"
    echo ""
    exit 1
fi

bold=$(remove_ext $1)

cp $FSLDIR/data/standard/MNI152_T1_1mm.nii.gz .
t1=MNI152_T1_1mm

thresh=10
conn=6

bet $t1 ${t1}_brain -R -m
fslmaths $t1 -subsamp2 ${t1}_2mm
fslmaths ${t1}_brain -subsamp2 ${t1}_brain_2mm

fslmaths $bold.nii.gz -Tmean func_mean

bold=func_mean

bet $bold ${bold}_brain -m -R

epi_reg --epi=${bold}_brain.nii.gz --t1=${t1}_2mm.nii.gz --t1brain=${t1}_brain_2mm.nii.gz --out=BOLD-to-T1 #registration
slices ${t1}_2mm.nii.gz BOLD-to-T1.nii.gz #make sure to check!

convert_xfm -omat T1-to-BOLD.mat -inverse BOLD-to-T1.mat #inverse matrix from bold-to-T1 to T1-to-bold
flirt -in ${t1}_brain_2mm.nii.gz -ref ${bold}.nii.gz -out T1-to-BOLD -init T1-to-BOLD.mat -applyxfm #register T1 to BOLD space

fast -g -o fast T1-to-BOLD.nii.gz #brain segmentation (white/grey/csf)

fslmaths fast_pve_0.nii.gz -thr 1 csf
fslmaths fast_pve_1.nii.gz -thr 1 greymatter
fslmaths fast_pve_2.nii.gz -thr 1 whitematter

fslmaths ${bold}_brain_mask.nii.gz -ero ${bold}_brain_mask_ero

fslmaths whitematter.nii.gz -mul ${bold}_brain_mask_ero.nii.gz whitematter
fslmaths greymatter.nii.gz -mul ${bold}_brain_mask_ero.nii.gz greymatter

cluster -i greymatter.nii.gz -t 1 --connectivity=6 --no_table --osize=greymatter_cluster_size
fslmaths greymatter_cluster_size.nii.gz -thr 10 -bin greymatter -odt char
imrm greymatter_cluster_size

cluster -i whitematter.nii.gz -t 1 --connectivity=6 --no_table --osize=whitematter_cluster_size
fslmaths whitematter_cluster_size.nii.gz -thr 10 -bin whitematter -odt char
imrm whitematter_cluster_size

rm BOLD-to* fast_* T1-to* ${bold}_brain* ${t1}_*
