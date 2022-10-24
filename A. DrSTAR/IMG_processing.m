%%% Function created by Nawara T. (Mattheyses lab - 05/20/2022) compatible
%%% with MATLAB R2020b / Function allows user to pic up background ROI for
%%% each cell, perfroms backgrond subtraction, flat field correctoin and
%%% bleach correction (simple ratio)

%This work is licensed under the Creative Commons Attribution 4.0
%International License. To view a copy of this license, visit
%http://creativecommons.org/licenses/by/4.0/ or send a letter to Creative
%Commons, PO Box 1866, Mountain View, CA 94042, USA.

function IMG_processing(yourpath, avg_FFC_488, avg_FFC_647)

    max_FFC_488 = max(avg_FFC_488(:));
    max_FFC_647 = max(avg_FFC_647(:));
    ContentInFold = dir([yourpath, '/Data_analysis']);
    background_488 = zeros(length(ContentInFold),1);
    background_647 = zeros(length(ContentInFold),1);
    
    %Pick background ROIs for a cellthen press OK 
    for i = 1:length(ContentInFold)
        if strfind(ContentInFold(i).name, 'STAR') >= 1
            T488_raw_path = [yourpath, '/', 'Data_analysis', '/', ContentInFold(i).name, '/', ContentInFold(i).name, '_488_Raw.tif'];
            T647_raw_path = [yourpath, '/', 'Data_analysis', '/', ContentInFold(i).name, '/', ContentInFold(i).name, '_647_Raw.tif'];
            r_488 = bfGetReader(T488_raw_path);
            I_488 = bfGetPlane(r_488, 1);
            r_647 = bfGetReader(T647_raw_path);
            I_647 = bfGetPlane(r_647, 1);
            
            close all
            figure(1)
            imagesc(I_488); axis equal
            roi1 = drawrectangle;
            roi1.Position = round(roi1.Position);
            choice = menu('Press OK after picking background','OK');  
            background_488(i,1) = mean(I_488(roi1.Position(2):(roi1.Position(2)+roi1.Position(4)),roi1.Position(1):(roi1.Position(1)+roi1.Position(3))), 'all');
            background_647(i,1) = mean(I_647(roi1.Position(2):(roi1.Position(2)+roi1.Position(4)),roi1.Position(1):(roi1.Position(1)+roi1.Position(3))), 'all'); 
            close all
        end
    end
fprintf('Background selection compleated')

    parfor (i = 1:length(ContentInFold))
        if strfind(ContentInFold(i).name, 'STAR') >= 1 
            %488 processing
            T488_raw_path = [yourpath, '/', 'Data_analysis', '/', ContentInFold(i).name, '/', ContentInFold(i).name, '_488_Raw.tif'];
            T488_raw = bfopen(T488_raw_path);
            T488_cor = cell(size(T488_raw{1,1},1),1);%fixed length2size 020322 will work 1 picture images   
            
            for ii = 1:size(T488_raw{1,1},1) %fixed length2size 020322 will work 1 picture images   
                T488_cor{ii,1} = zeros(size(T488_raw{1,1}{1,1}));
            end
                
            for ii = 1:size(T488_raw{1,1},1) %fixed length2size 020322 will work 1 picture images   
                
                T488_raw{1,1}{ii, 1} = T488_raw{1,1}{ii, 1} - background_488(i,1); %Background subtraction
                mask_488 = T488_raw{1,1}{ii, 1} >= 0;
                T488_raw{1,1}{ii, 1} = double(T488_raw{1,1}{ii, 1}) .* mask_488; %preserving non negative values only as in int16-bit image
                T488_raw{1,1}{ii, 1} = (double(T488_raw{1,1}{ii, 1}) ./ double(avg_FFC_488)).* double(max_FFC_488); %Flat field correction
                ratio = mean(T488_raw{1,1}{1, 1}, 'all') / mean(T488_raw{1,1}{ii, 1}, 'all');
                T488_cor{ii, 1} = T488_raw{1,1}{ii, 1} .* ratio; %Bleach correction
            end
                 
            fn488 = [yourpath, '/', 'Data_analysis', '/', ContentInFold(i).name, '/', ContentInFold(i).name, '_488_Cor.tif'];
            cell_mat2tiff(fn488, T488_cor)
            T488_raw = [];
            T488_cor = [];
            
        
            %647 processing
            T647_raw_path = [yourpath, '/', 'Data_analysis', '/', ContentInFold(i).name, '/', ContentInFold(i).name, '_647_Raw.tif'];
            T647_raw = bfopen(T647_raw_path);
            T647_cor = cell(size(T647_raw{1,1},1),1); %fixed length2size 020322 will work 1 picture images   
            
            for ii = 1:size(T647_raw{1,1},1)
                T647_cor{ii,1} = zeros(size(T647_raw{1,1}{1,1})); %fixed length2size 020322 will work 1 picture images   
            end
            
            for ii = 1:size(T647_raw{1,1},1) %fixed length2size 020322 will work 1 picture images   
                
                T647_raw{1,1}{ii, 1} = T647_raw{1,1}{ii, 1} - background_647(i,1); %Background subtraction
                mask_647 = T647_raw{1,1}{ii, 1} >= 0;
                T647_raw{1,1}{ii, 1} = double(T647_raw{1,1}{ii, 1}) .* mask_647; %preserving non negative values only as in int16-bit image
                T647_raw{1,1}{ii, 1} = (double(T647_raw{1,1}{ii, 1}) ./ double(avg_FFC_647)).* double(max_FFC_647); %Flat field correction
                ratio = mean(T647_raw{1,1}{1, 1}, 'all') / mean(T647_raw{1,1}{ii, 1}, 'all');
                T647_cor{ii, 1} = T647_raw{1,1}{ii, 1} .* ratio; %Bleach correction
            end
            
            fn647 = [yourpath, '/', 'Data_analysis', '/', ContentInFold(i).name, '/', ContentInFold(i).name, '_647_Cor.tif'];
            cell_mat2tiff(fn647, T647_cor)
            T647_raw = [];
            T647_cor = [];    
        end
    end
    writematrix([], [yourpath, '/', 'Corrections', '/','Corrections_done.txt']); %added 032122 for automation if somehting gees wrong
end