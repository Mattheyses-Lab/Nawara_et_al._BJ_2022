%%% Function created by Nawara T. (Mattheyses lab - 05/20/2022) compatible
%%% with MATLAB R2020b / Function uses the generated caliblration file and
%%% performs AVG intesity porjection for FFC files splits them based on
%%% emisison and applys a blur to them defined by sigma

%This work is licensed under the Creative Commons Attribution 4.0
%International License. To view a copy of this license, visit
%http://creativecommons.org/licenses/by/4.0/ or send a letter to Creative
%Commons, PO Box 1866, Mountain View, CA 94042, USA.

function [avg_FFC_488, avg_FFC_647] = FFC_splitter(yourpath, OS_Calibration_File)

ContentInFold = dir(yourpath);
sigma = 10;

% Loop on each folder
    for i = 1:length(ContentInFold) % start at 3 to skip . and .. and DS_store
        if strfind(ContentInFold(i).name, 'FF488') >= 1     
            image_to_be_splitted = fullfile(yourpath,ContentInFold(i).name);
            whole_img = bfopen(image_to_be_splitted);
            %prealocation of matrixes for splitted image
            T488_raw = cell(length(whole_img{1, 1}),1);
            
            %FFC 488 processing
            for ii = 1:length(whole_img{1, 1}) %splitting image based on ROI data from OS_cali.m
                T488_raw{ii,1} = uint16(zeros(OS_Calibration_File(1,4)+1, OS_Calibration_File(1,3)+1));
                T488_raw{ii,1} = whole_img{1,1}{ii,1}(OS_Calibration_File(1,2):(OS_Calibration_File(1,2)+OS_Calibration_File(1,4)),OS_Calibration_File(1,1):(OS_Calibration_File(1,1)+OS_Calibration_File(1,3)));
            end
            
            T488_raw = cat(3,T488_raw{:});
            avg_FFC_488 = mean(T488_raw,3);
            avg_FFC_488 = imgaussfilt(avg_FFC_488, sigma);  %added as the uSlides have some pores 020822 TN
            
            fn488 = [yourpath, '/', 'Corrections', '/', ContentInFold(i).name(1:end-8), '.tif']; %032122 chaged the saving name TN
            cell_mat2tiff(fn488, avg_FFC_488)
            %save(fn488,'avg_FFC_488', '-v7.3');
            
        elseif strfind(ContentInFold(i).name, 'FF647') >= 1
            image_to_be_splitted = fullfile(yourpath,ContentInFold(i).name);
            whole_img = bfopen(image_to_be_splitted);
            %prealocation of matrixes for splitted image
            T647_raw = cell(length(whole_img{1, 1}),1);
           
            %FFC 647 processing
            for ii = 1:length(whole_img{1, 1}) %splitting image based on ROI data from OS_cali.m
                T647_raw{ii,1} = zeros(OS_Calibration_File(2,4)+1, OS_Calibration_File(2,3)+1);
                T647_raw{ii,1} = whole_img{1,1}{ii,1}(OS_Calibration_File(2,2):(OS_Calibration_File(2,2)+OS_Calibration_File(2,4)),OS_Calibration_File(2,1):(OS_Calibration_File(2,1)+OS_Calibration_File(2,3)));              
            end
            
            T647_raw = cat(3,T647_raw{:});
            avg_FFC_647 = mean(T647_raw,3);
            avg_FFC_647 = imgaussfilt(avg_FFC_647, sigma); %added as the uSlides have some pores 020822 TN
            
            fn647 = [yourpath, '/', 'Corrections', '/', ContentInFold(i).name(1:end-8), '.tif']; %032122 chaged the saving name TN
            cell_mat2tiff(fn647, avg_FFC_647)
            %save(fn647,'avg_FFC_647', '-v7.3');
        end    
    end  
end