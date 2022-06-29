%This work is licensed under the Creative Commons Attribution 4.0
%International License. To view a copy of this license, visit
%http://creativecommons.org/licenses/by/4.0/ or send a letter to Creative
%Commons, PO Box 1866, Mountain View, CA 94042, USA.

function IMGREG_data(yourpath)
%%% Created by Nawara T. Matttheyses Lab @UAB 03.02.20
%IMREG is a fucntion that allows to register multichannel image stack,
%based on the calibration grid. it needs 7 types of imputs:
% 1. Path to the referance grid picture [path488];
% 2. Path to to be registered grid picture [path647];
% 3. Path to the image or image stack to be registered [pathdata_reg];
% 4. Path to the referance image or image stack [pathdata_ref];
% 5. Specifed sigma for Gausian filtering
% 6. Defined 'FilterSize'/f for feature detection
% 7. Name a of the file to be saved
% Program requires instalation of bio-formats tool box details:
% https://docs.openmicroscopy.org/bio-formats/6.4.0/users/matlab/index.html
% Requires Computer Vision Toolbox
% For any qestions email tomasz.nawara.95@gmail.com

%[last edited:05.19.22]
%  Fixed error in how detected points were saved in matrix [04.26.20]
%  Fixed image display [04.27.20]
%  Added the registration conformation window (figure 1) [05.05.20]
%  Added Registration of signle frame images [07.04.20]
%  Added batch processing [10.09.20]
%  Added batch processing for single images [11.10.20]
%  Display for registration confomration [04.11.22]
%  Changed registration for loop to 3D matrix [05.19.22]

ContentInFold = dir(yourpath);

for i = 1:length(ContentInFold)
    if strfind(ContentInFold(i).name, 'STAR') >= 1
        T647_cor_path = [yourpath, '/', 'Data_analysis', '/', ContentInFold(i).name(1:end-4), '/', ContentInFold(i).name(1:end-4), '_647_Cor.tif'];
        fn647 = [yourpath, '/', 'Data_analysis', '/', ContentInFold(i).name(1:end-4), '/', ContentInFold(i).name(1:end-4), '_647_Cor_Reg.tif']; %%%%
        
        if ~exist(fn647, 'file') %%%%
                 
            T647_Cor = bfopen(T647_cor_path);
            T647_Cor_Reg = cell(size(T647_Cor{1,1},1),1);  %fixed length2size 020322 will work 1 picture images
            
            T488_cor_path = [yourpath, '/', 'Data_analysis', '/', ContentInFold(i).name(1:end-4), '/', ContentInFold(i).name(1:end-4), '_488_Cor.tif'];
            r_488 = bfGetReader(T488_cor_path);
            I_488 = bfGetPlane(r_488, 1);
            
            for ii = 1:size(T647_Cor{1,1},1)
                T647_Cor_Reg{ii,1} = zeros(size(T647_Cor{1,1}{1,1}));
            end
            
            [Ireg,O_trans,Spacing,~,~,~] = image_registration(double(T647_Cor{1,1}{1,1}),double(I_488));
            
            %figure; imshowpair(I_488*200, T647_Cor{1,1}{1,1}*200);
            %figure; imshowpair(I_488*200, uint16(Ireg)*200);
            figure; imshow(double(T647_Cor{1,1}{1,1}) ./ imgaussfilt(double(I_488),1), [0 5])
            
            for k = 1:size(T647_Cor_Reg,1)
                 T647_Cor_Reg{k, 1} = bspline_transform(O_trans,double(T647_Cor{1,1}{k,1}),Spacing,2);
            end
            
            cell_mat2tiff(fn647, T647_Cor_Reg);
            T647_Cor_Cor = [];
            T647_Cor_Reg = [];
        end %%$%%
    end
end
writematrix([], [yourpath, '/', 'Corrections', '/','647_Registration_done.txt']); %added 032122 for automation if somehting gees wrong
fprintf('<3 Registration completed <3 \n');
end

