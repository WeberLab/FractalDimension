csfmask = niftiread('csf_mask.nii.gz');
csflog = logical(csfmask);
greymask = niftiread('grey_mask.nii.gz');
greylog = logical(greymask);
whitemask = niftiread('white_mask.nii.gz');
whitelog = logical(whitemask);
noisemask = niftiread('noise_mask.nii.gz');
noiselog = logical(noisemask);
sub8ins = niftiread('sub-08_ses-1_task-ins_bold.nii.gz');



[row,col,page] = ind2sub(size(greylog),datasample(find(greylog),1))
samplegrey = squeeze(sub8ins(row,col,page,:));


[row,col,page] = ind2sub(size(noiselog),datasample(find(noiselog),1))
samplenoise = squeeze(sub8ins(row,col,page,:));
HurstOlga(samplenoise(3:end),2.5);
