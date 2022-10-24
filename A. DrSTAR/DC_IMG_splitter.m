%%% Function created by Nawara T. (Mattheyses lab - 05/20/2022) compatible
%%% with MATLAB R2020b / This is like Cairn Image Splitter - splits
%%% live-cell data based on emision wavelght defined by OS_Calibration file

%This work is licensed under the Creative Commons Attribution 4.0
%International License. To view a copy of this license, visit
%http://creativecommons.org/licenses/by/4.0/ or send a letter to Creative
%Commons, PO Box 1866, Mountain View, CA 94042, USA.

function DC_IMG_splitter(yourpath, Calibration_File)

ContentInFold = dir(yourpath);
    parfor i = 1:length(ContentInFold) % start at 3 to skip . and .. and DS_store
        if strfind(ContentInFold(i).name, 'STAR_488') >= 1
            dir4cell = [yourpath, '/', 'Data_analysis', '/', ContentInFold(i).name(1:end-4)];
            index4save = strfind(dir4cell,'488');
            mkdir([dir4cell(1:index4save-1) dir4cell(index4save+4:end)]); warning('off');
            image_to_be_splitted_488 = fullfile(yourpath,ContentInFold(i).name);
            index4lambda = strfind(image_to_be_splitted_488,'488');
            image_to_be_splitted_647 = image_to_be_splitted_488;
            image_to_be_splitted_647(index4lambda:index4lambda+2) = '647';
     
                whole_img_488 = bfopen(image_to_be_splitted_488);
                whole_img_647 = bfopen(image_to_be_splitted_647);
                %prealocation of matrixes for splitted image
                T488_raw = cell(size(whole_img_488{1, 1},1),1);  %fixed length2size 020322 will work 1 picture images   
                T647_raw = cell(size(whole_img_647{1, 1},1),1);  %fixed length2size 020322 will work 1 picture images 

                for ii = 1:size(whole_img_488{1, 1},1) %fixed length2size 020322 will work 1 picture images  
                    T488_raw{ii,1} = uint16(zeros(Calibration_File(1,4)+1, Calibration_File(1,3)+1));
                end
                
                for ii = 1:size(whole_img_488{1, 1},1) %splitting image based on ROI data from OS_cali.m %fixed length2size 020322 will work 1 picture images  
                    T488_raw{ii,1} = whole_img_488{1,1}{ii,1}(Calibration_File(1,2):(Calibration_File(1,2)+Calibration_File(1,4)),Calibration_File(1,1):(Calibration_File(1,1)+Calibration_File(1,3))); 
                end
                
                cell_name = ContentInFold(i).name(1:end-4);
                index4cell_name = strfind(cell_name,'488');
                cell_name = [cell_name(1:index4cell_name-1) cell_name(index4cell_name+4:end)];
                
                fn488 = [dir4cell(1:index4save-1) dir4cell(index4save+4:end), '/', cell_name, '_488_Raw.tif'];
                cell_mat2tiff(fn488, T488_raw)
                T488_raw = [];
                
                for ii = 1:size(whole_img_647{1, 1},1) %fixed length2size 020322 will work 1 picture images  
                    T647_raw{ii,1} = uint16(zeros(Calibration_File(2,4)+1, Calibration_File(2,3)+1));
                end
                
                for ii = 1:size(whole_img_647{1, 1},1) %splitting image based on ROI data from OS_cali.m %fixed length2size 020322 will work 1 picture images  
                    T647_raw{ii,1} = whole_img_647{1,1}{ii,1}(Calibration_File(2,2):(Calibration_File(2,2)+Calibration_File(2,4)),Calibration_File(2,1):(Calibration_File(2,1)+Calibration_File(2,3))); 
                end

                fn647 = [dir4cell(1:index4save-1) dir4cell(index4save+4:end), '/', cell_name, '_647_Raw.tif'];
                cell_mat2tiff(fn647, T647_raw)
                T647_raw = [];
    %             save(fn488,'T488_raw', '-v7.3');
    %             save(fn647,'T488_raw', '-v7.3'); 
        end
    end
    writematrix([], [yourpath, '/', 'Corrections', '/','Splitting_done.txt']); %added 032122 for automation if somehting gees wrong
end

