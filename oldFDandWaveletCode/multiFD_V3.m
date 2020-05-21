% Once excel data is cut to include only exams that we want FD for and
% those exams have appropriate MRS X Y data, open master_table_maker_V2 and
% change the excel file name to the appropriate file


master_table_maker_V2   % use master_table_maker_V2 which does not extract FD data
cd FD

for a=1:size(masterdata.exams_done.data,1)
    
    if isnan (masterdata.exams_done.data(a,5))
        
        % Do nothing
        
    else
        
        examID=(masterdata.exams_done.data(a,1));
        
        MRSx=(masterdata.MRS_voxel.data(a,2));
        
        MRSy=(masterdata.MRS_voxel.data(a,3));
        
        FDall2_V5_MRSroi
        
    end
    
end

        
cd ..        
        