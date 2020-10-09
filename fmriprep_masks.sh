#!/bin/bash

# Author: Alex Weber, Ph.D.
# Date: Oct 9 2020

if [ $# -lt 1 ]
then
    echo "Usage: $0 <fmriprepBOLD brain mask>"
    echo "Example: $0 sub-01_task-rest_space-MNI152NLin2009cAsym_desc-brain_mask.nii.gz"
    echo "This grabs the T1 from MNI space, then creates very conservative whitematter and greymatter masks"
    echo "[note: default cluster removal is 10 and 6]"
    echo ""
    exit 1
fi

boldmask=$(remove_ext $1)

cp $FSLDIR/data/standard/MNI152_T1_1mm.nii.gz .

t1=MNI152_T1_1mm.nii.gz

bet $t1 ${t1}_brain -R -m

fslmaths ${t1}_brain -subsamp2 ${t1}_brain_2mm

fast -g -o fast ${t1}_brain_2mm

fslmaths fast_pve_0.nii.gz -thr 1 csf
fslmaths fast_pve_1.nii.gz -thr 1 greymatter
fslmaths fast_pve_2.nii.gz -thr 1 whitematter

fslmaths ${boldmask} -ero ${boldmask}_ero

fslmaths whitematter.nii.gz -mul ${boldmask}_ero.nii.gz whitematter
fslmaths greymatter.nii.gz -mul ${boldmask}_ero.nii.gz greymatter

cluster -i greymatter.nii.gz -t 1 --connectivity=6 --no_table --osize=greymatter_cluster_size
fslmaths greymatter_cluster_size.nii.gz -thr 10 -bin greymatter -odt char
imrm greymatter_cluster_size

cluster -i whitematter.nii.gz -t 1 --connectivity=6 --no_table --osize=whitematter_cluster_size
fslmaths whitematter_cluster_size.nii.gz -thr 10 -bin whitematter -odt char
imrm whitematter_cluster_size

rm fast_* ${boldmask}_ero* ${t1}_*
