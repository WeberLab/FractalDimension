%This script loads a nifti sample file, applies a grey mask, and calls FracTool function

%load nifti file and mask
sample_file = niftiread('../sub-08_ses-1/func/sub-08_ses-1_task-ins_bold.nii.gz');
grey_mask = niftiread('../sub-08_ses-1/func/grey_mask.nii.gz');
grey_log = logical(grey_mask);
[row,col,page] = ind2sub(size(grey_log),datasample(find(grey_log),1))
%get data from the random row, col, page
sample_grey = squeeze(sample_file(row,col,page,:));
%convert to type double
test_signal = double(sample_grey);
%function call 
output = Fractool(test_signal)
%a file 'FracTool.txt' is generated in the directory
