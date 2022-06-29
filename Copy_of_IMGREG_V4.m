%This work is licensed under the Creative Commons Attribution 4.0
%International License. To view a copy of this license, visit
%http://creativecommons.org/licenses/by/4.0/ or send a letter to Creative
%Commons, PO Box 1866, Mountain View, CA 94042, USA.

function IMGREG_V4(yourpath, Grid_488_raw, Grid_647_reg, use_fixed_detection)
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

sigma = 3;
f = 11;
zz = 0;
what_a_day_to_code = 0;
dbstop in IMGREG_V2 at 131 if m=='N'

ContentInFold = dir(yourpath);

for i = 1:length(ContentInFold)
    if strfind(ContentInFold(i).name, 'STAR') >= 1 & what_a_day_to_code == 0
        T488_cor_path = [yourpath, '/', 'Data_analysis', '/', ContentInFold(i).name(1:end-4), '/', ContentInFold(i).name(1:end-4), '_488_Cor.tif'];
        T647_cor_path = [yourpath, '/', 'Data_analysis', '/', ContentInFold(i).name(1:end-4), '/', ContentInFold(i).name(1:end-4), '_647_Cor.tif'];
        fn647 = [yourpath, '/', 'Data_analysis', '/', ContentInFold(i).name(1:end-4), '/', ContentInFold(i).name(1:end-4), '_647_Cor_Reg.tif'];
        tform_path = [yourpath, '/', 'Data_analysis', '/', ContentInFold(i).name(1:end-4), '/' ContentInFold(i).name(1:end-4), '_tform.mat'];
        fig_path = [yourpath, '/', 'Data_analysis', '/', ContentInFold(i).name(1:end-4), '/' ContentInFold(i).name(1:end-4), '_reg_fig.png'];
        
        if exist(tform_path ,'file')
            load(tform_path); %%%%
            I = size(Grid_488_raw);
            
        else
            I = size(Grid_488_raw);
            r_488 = bfGetReader(T488_cor_path);
            I_488 = bfGetPlane(r_488, 1);
            r_647 = bfGetReader(T647_cor_path);
            I_647 = bfGetPlane(r_647, 1);
            I_647_raw = I_647;
            
            [tform, p4p_488, p4p_647] = tform_generation(I_488, I_647);
            I_647 = imwarp(I_647, tform, 'OutputView', imref2d(I));
            
            figure('units','normalized','outerposition',[0 0 1 1])
            tiledlayout(2,2);
            nexttile; imshowpair(imadjust(I_488), imadjust(I_647_raw), 'ColorChannels', 'green-magenta'); title('Unregistered overlay'); hold on
            plot(p4p_488(1:2:end),p4p_488(2:2:end),'g+', 'MarkerSize', 10)
            plot(p4p_647(1:2:end),p4p_647(2:2:end),'r+', 'MarkerSize', 10)
            nexttile; imshowpair(imadjust(I_488), imadjust(I_647), 'ColorChannels', 'green-magenta'); title('Registered overlay'); hold on
            nexttile; imshow(double(I_647_raw)./double(I_488), [0 5]); title('Raw 647/488 ratio'); hold on; set(gca,'xticklabel',{[]}); set(gca,'yticklabel',{[]})
            nexttile; imshow(double(I_647)./double(imgaussfilt(I_488,1)), [0 5]); title('Reg 647/488 ratio'); hold on; set(gca,'xticklabel',{[]}); set(gca,'yticklabel',{[]})
            
            m = input('Do you want to continue, Y/N:','s');
            
            if m=='Y'
                
                FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
                save_tfrom(tform_path, tform, fig_path, FigList)
                
            end
        end
        
        T647_Cor = bfopen(T647_cor_path);
        T647_Cor_Reg = cell(size(T647_Cor{1,1},1),1);  %fixed length2size 020322 will work 1 picture images
        
        for ii = 1:size(T647_Cor{1,1},1)
            T647_Cor_Reg{ii,1} = zeros(size(T647_Cor{1,1}{1,1}));
        end
        
        %%% Depracted
        %             for k = 1:size(T647_Cor_Reg,1)
        %                 T647_Cor_Reg{k, 1} = imwarp(T647_Cor{1,1}{k,1}, tform, 'OutputView', imref2d(I));
        %             end
        
        AA = uint16(zeros([size(T647_Cor{1,1}{1,1}), size(T647_Cor{1,1},1)]));
        BB = AA;
        for dd = 1:size(T647_Cor{1,1},1)
            AA(:,:,dd) = T647_Cor{1,1}{dd,1};
        end
        
        BB(:,:,:) = imwarp(AA(:,:,:), tform, 'OutputView', imref2d(I)); %registration
        
        for dd = 1:size(T647_Cor{1,1},1)
            T647_Cor_Reg{dd, 1} = BB(:,:,dd);
        end
        
        %fn647 = [yourpath, '/', 'Data_analysis', '/',
        %ContentInFold(i).name(1:end-4), '/',
        %ContentInFold(i).name(1:end-4), '_647_Cor_Reg.tif']; %%%%
        cell_mat2tiff(fn647, T647_Cor_Reg);
        T647_Cor_Cor = [];
        T647_Cor_Reg = [];
    end
end
end
%
%
%                 %%%%%%%%%%%%%%
%
%
%
%                 if use_fixed_detection == 0
%                     %Apply Gausian smotching to pictures
%                     I488 = imgaussfilt(Grid_488_raw, sigma);
%                     I647 = imgaussfilt(Grid_647_reg, sigma);
%
%                     %Use bfopen to load the data
%                     I = size(Grid_488_raw); %size of the referance picture
%
%                     %Detect features
%                     points488_raw = detectHarrisFeatures(I488, 'FilterSize', f);
%                     points488 = sortrows(points488_raw.Location, [-2 1]);
%
%                     points647_raw = detectHarrisFeatures(I647, 'FilterSize', f);
%                     points647 = sortrows(points647_raw.Location, [-2 1]);
%
%                     %     figure(1)
%                     %     tiledlayout(2,2)
%                     %     nexttile; imagesc(Grid_488_raw); title('Raw 488')
%                     %     nexttile; imshow(imadjust(I488)); hold on; plot(points488_raw); title('Detection 488'); hold off;
%                     %     nexttile; imagesc(Grid_647_reg); title('Raw 647')
%                     %     nexttile; imshow(imadjust(I488)); hold on; plot(points488_raw); title('Detection 647'); hold off;
%
%                     sorted_488 = zeros(length(points488), 2);
%                     sorted_647 = zeros(length(points647), 2);
%
%                     %Accept detection
%                     %     m=input('Do you want to continue, Y/N:','s');
%                     %     if m=='N'
%                     %         close all
%                     %         return
%                     %
%                     %     elseif m=='Y'
%                     %Create registration mask
%                     %but correct for terrible matching of
%                     %detection grid
%                     %if you delet numbers bigger from around xzero then you automatically
%                     %elimnaets unaligned points :) done
%                     for gg = 1:length(points488)
%                         if points488(gg,2) - points488(gg+1,2) > 8
%                             break
%                         end
%                     end
%
%                     for ff = 0:(length(points488)/gg)-1
%                         sorted_488(1+(gg*ff):gg+(gg*ff),:) = sortrows(points488(1+(gg*ff):gg+(gg*ff), :), -1);
%                         sorted_647(1+(gg*ff):gg+(gg*ff),:) = sortrows(points647(1+(gg*ff):gg+(gg*ff), :), -1);
%                     end
%
%                     n = length(sorted_488); % number of grid holes to be used for registration (x,y coordinates)
%                     tform = fitgeotrans(sorted_647, sorted_488, 'lwm', n);
%                     %     end
%                     close all
%
%                     ContentInFold = dir(yourpath);
%
%
%
%                     for i = 1:length(ContentInFold)
%                         if strfind(ContentInFold(i).name, 'STAR') >= 1 & what_a_day_to_code == 0
%                             T488_cor_path = [yourpath, '/', 'Data_analysis', '/', ContentInFold(i).name(1:end-4), '/', ContentInFold(i).name(1:end-4), '_488_Cor.tif'];
%                             T647_cor_path = [yourpath, '/', 'Data_analysis', '/', ContentInFold(i).name(1:end-4), '/', ContentInFold(i).name(1:end-4), '_647_Cor.tif'];
%                             r_488 = bfGetReader(T488_cor_path);
%                             I_488 = bfGetPlane(r_488, 1);
%                             r_647 = bfGetReader(T647_cor_path);
%                             I_647 = bfGetPlane(r_647, 1);
%                             I_647_raw = I_647;
%                             I_647 = imwarp(I_647, tform, 'OutputView', imref2d(I));
%
%                             if zz == 0
%
%                                 zz = 1;
%                                 figure(1);
%                                 tiledlayout(2,2)
%                                 nexttile; imshow(imadjust(Grid_488_raw)); title('Detection 488'); hold on; set(gca,'xticklabel',{[]}); set(gca,'yticklabel',{[]})
%                                 plot(points488(:,1),points488(:,2),'r+'); hold off;
%                                 nexttile; imshowpair(imadjust(I_488), imadjust(I_647_raw), 'ColorChannels', 'green-magenta'); title('Unregistered overlay');
%                                 nexttile; imshow(imadjust(Grid_647_reg)); title('Detection 647'); hold on; set(gca,'xticklabel',{[]}); set(gca,'yticklabel',{[]})
%                                 plot(points647(:,1),points647(:,2),'r+'); hold off;
%                                 nexttile; imshowpair(imadjust(I_488), imadjust(I_647), 'ColorChannels', 'green-magenta'); title('Registered overlay');
%                                 figure(2); imshow(double(I_647)./double(I_488), [0 5]); title('647/488 ratio'); hold on; set(gca,'xticklabel',{[]}); set(gca,'yticklabel',{[]})
%                                 m=input('Do you want to continue, Y/N:','s');
%
%                                 if m=='N'
%                                     close all; points488 = sorted_488; points647 = sorted_647;
%                                     openvar('points488');
%                                     openvar('points647'); diff488_647 = points488 - points647; display(diff488_647); openvar('diff488_647');
%                                     fprintf('Fix the detection points manually\n first examin the diff488_647 matrix\n if values are bigger than around 0\n open the point488 and point 647\n and match them by hand\n Once matched press continue\n')
%                                     n = length(points488); % number of grid holes to be used for registration (x,y coordinates)
%                                     %                     save([yourpath, '/Corrections/points488_fixed.mat'], 'points488')
%                                     %                     save([yourpath, '/Corrections/points647_fixed.mat'], 'points647')
%                                     tform = fitgeotrans(points647, points488, 'lwm', n);
%                                     I_647 = imwarp(I_647_raw, tform, 'OutputView', imref2d(I));
%                                     imshowpair(imadjust(I_488), imadjust(I_647), 'ColorChannels', 'green-magenta'); title('NEW Registered overlay'); axis equal;
%
%                                     m=input('Do you want to continue, Y/N:','s');
%                                     if m=='N'
%                                         fprintf('Translate 647 to match 488 perfetlly\n')
%                                         assignin('base','I_488_2bTRANS', I_488)
%                                         assignin('base','I_647_2bTRANS', I_647)
%                                         OS_IMG_TRANSLATION_V2
%                                         %open the calibration appp then extract X and Y
%                                         %and add it to the points 647 (translation)
%
%                                         choice = menu('Press OK after aligment','OK');
%
%                                         X_TRANS = evalin('base','X_TRANS');
%                                         Y_TRANS = evalin('base','Y_TRANS');
%
%                                         if isequal(X_TRANS, [])
%                                             X_TRANS = 0;
%                                         end
%
%                                         if isequal(Y_TRANS, [])
%                                             Y_TRANS = 0;
%                                         end
%
%                                         points647(:, 1) = points647(:, 1) -  X_TRANS;
%                                         points647(:, 2) = points647(:, 2) -  Y_TRANS;
%
%                                         tform = fitgeotrans(points647, points488, 'lwm', n);
%                                         I_647 = imwarp(I_647_raw, tform, 'OutputView', imref2d(I));
%                                         imshowpair(imadjust(I_488), imadjust(I_647), 'ColorChannels', 'green-magenta'); title('NEW Registered overlay'); axis equal;
%
%                                         m=input('Do you want to continue, Y/N:','s');
%                                         if m=='N'
%                                             fprintf('Registration mask generation failed')
%                                             return
%                                         end
%
%                                         save([yourpath, '/Corrections/points488_fixed.mat'], 'points488')
%                                         save([yourpath, '/Corrections/points647_fixed.mat'], 'points647')
%
%                                     elseif m=='Y'
%                                         save([yourpath, '/Corrections/points488_fixed.mat'], 'points488')
%                                         save([yourpath, '/Corrections/points647_fixed.mat'], 'points647')
%
%                                     end
%
%
%                                 elseif m=='Y'
%                                     close all
%                                 end
%
%                             end
%                             what_a_day_to_code = 1;
%                             save([yourpath, '/Corrections/tform.mat'], 'tform'); %%%%
%                         end
%                     end
%                 else %after points fixetaion just load inthe fixed points
%                     ContentInFold = dir(yourpath);
%                     I = size(Grid_488_raw);
%                     load([yourpath, '/Corrections/points488_fixed.mat'])
%                     load([yourpath, '/Corrections/points647_fixed.mat'])
%                     n = length(points488);
%                     tform = fitgeotrans(points647, points488, 'lwm', n);
%                     save([yourpath, '/Corrections/tform.mat'], 'tform'); %%%%
%                 end
%             end %%%%
%
%             parfor i = 1:length(ContentInFold)
%                 if strfind(ContentInFold(i).name, 'STAR') >= 1
%                     T647_cor_path = [yourpath, '/', 'Data_analysis', '/', ContentInFold(i).name(1:end-4), '/', ContentInFold(i).name(1:end-4), '_647_Cor.tif'];
%                     fn647 = [yourpath, '/', 'Data_analysis', '/', ContentInFold(i).name(1:end-4), '/', ContentInFold(i).name(1:end-4), '_647_Cor_Reg.tif']; %%%%
%
%                     if ~exist(fn647, 'file') %%%%
%
%                         T647_Cor = bfopen(T647_cor_path);
%                         T647_Cor_Reg = cell(size(T647_Cor{1,1},1),1);  %fixed length2size 020322 will work 1 picture images
%
%                         for ii = 1:size(T647_Cor{1,1},1)
%                             T647_Cor_Reg{ii,1} = zeros(size(T647_Cor{1,1}{1,1}));
%                         end
%
%                         %%% Depracted
%                         %             for k = 1:size(T647_Cor_Reg,1)
%                         %                 T647_Cor_Reg{k, 1} = imwarp(T647_Cor{1,1}{k,1}, tform, 'OutputView', imref2d(I));
%                         %             end
%
%                         AA = uint16(zeros([size(T647_Cor{1,1}{1,1}), size(T647_Cor{1,1},1)]));
%                         BB = AA;
%                         for dd = 1:size(T647_Cor{1,1},1)
%                             AA(:,:,dd) = T647_Cor{1,1}{dd,1};
%                         end
%
%                         BB(:,:,:) = imwarp(AA(:,:,:), tform, 'OutputView', imref2d(I)); %registration
%
%                         for dd = 1:size(T647_Cor{1,1},1)
%                             T647_Cor_Reg{dd, 1} = BB(:,:,dd);
%                         end
%
%                         %fn647 = [yourpath, '/', 'Data_analysis', '/',
%                         %ContentInFold(i).name(1:end-4), '/',
%                         %ContentInFold(i).name(1:end-4), '_647_Cor_Reg.tif']; %%%%
%                         cell_mat2tiff(fn647, T647_Cor_Reg);
%                         T647_Cor_Cor = [];
%                         T647_Cor_Reg = [];
%                     end %%$%%
%                 end
%             end
%             writematrix([], [yourpath, '/', 'Corrections', '/','647_Registration_done.txt']); %added 032122 for automation if somehting gees wrong
%             fprintf('<3 Registration completed <3 \n');
%         end
%     end
% end
% end

function save_tfrom(tform_path, tform, fig_path, FigList)
save(tform_path, 'tform');
saveas(FigList, fig_path)
end

function [tform, p_488, p_647] = tform_generation(I_488, I_647)

[Gmag,~] = imgradient(I_488*5000);
G = Gmag > 10;
BW_cell = imfill(~G,'holes');  %figure; imshow(BW_cell)
BW_cell = imerode(BW_cell, ones(10,10));
BW_cell_exlud = bwpropfilt(BW_cell,'Area',[20 10*10^7]); %figure; imshow(BW_G_exlud)
Cell_mask = imdilate(BW_cell_exlud, ones(10,10));

Gaus_1_488 = imgaussfilt(I_488, 2, 'FilterSize',7); %figure; imagesc(Gaus_1, [0 1000])
Gaus_2_488 = imgaussfilt(I_488, 4, 'FilterSize',7); %figure; imagesc(Gaus_2, [0 1000])
Gaus_dif_488 = (Gaus_1_488 - Gaus_2_488) .* uint16(Cell_mask); %figure; imagesc(Gaus_dif, [0 5])
Gaus_dif_smooth_488 = imgaussfilt(Gaus_dif_488,1); %figure; imagesc(Gaus_dif_smooth, [0 5])

Gaus_dif_NaN_488 = double(Gaus_dif_smooth_488);
Gaus_dif_NaN_488(Gaus_dif_NaN_488 < 0) = 0;
Gaus_dif_NaN_488(Gaus_dif_NaN_488 == 0) = NaN;

vec_Gaus_dif_NaN_488 = Gaus_dif_NaN_488;
vec_Gaus_dif_488 = vec_Gaus_dif_NaN_488(~isnan(vec_Gaus_dif_NaN_488));

mode_Gaus_dif_488 = mode(vec_Gaus_dif_488);

%Mask generation
BW_488 = Gaus_dif_NaN_488 > mode_Gaus_dif_488 ; %figure; imshow(BW);
BW2_488_size = bwpropfilt(BW_488, 'Area', [20 200]); %figure; imshow(BW2_488_size);
BW3_488 = bwpropfilt(BW2_488_size, "Eccentricity", [0 0.55]); %figure; imshow(BW3_488);
%BW3_488 = imdilate(BW2_488, ones(3,3)); figure; imshow(BW3_488);

Gaus_1_647 = imgaussfilt(I_647, 2, 'FilterSize',7); %figure; imagesc(Gaus_1, [0 1000])
Gaus_2_647 = imgaussfilt(I_647, 4, 'FilterSize',7); %figure; imagesc(Gaus_2, [0 1000])
Gaus_dif_647 = (Gaus_1_647 - Gaus_2_647) .* uint16(Cell_mask); %figure; imagesc(Gaus_dif, [0 5])
Gaus_dif_smooth_647 = imgaussfilt(Gaus_dif_647,1); %figure; imagesc(Gaus_dif_smooth, [0 5])

Gaus_dif_NaN_647 = double(Gaus_dif_smooth_647);
Gaus_dif_NaN_647(Gaus_dif_NaN_647 < 0) = 0;
Gaus_dif_NaN_647(Gaus_dif_NaN_647 == 0) = NaN;

vec_Gaus_dif_NaN_647 = Gaus_dif_NaN_647;
vec_Gaus_dif_647 = vec_Gaus_dif_NaN_647(~isnan(vec_Gaus_dif_NaN_647));

mode_Gaus_dif_647 = mode(vec_Gaus_dif_647);

%Mask generation
BW_647 = Gaus_dif_NaN_647 > mode_Gaus_dif_647 ; %figure; imshow(BW_647);
BW2_647_size = bwpropfilt(BW_647, 'Area', [20 200]); %figure; imshow(BW2_647_size);
BW3_647 = bwpropfilt(BW2_647_size, "Eccentricity", [0 0.55]); %figure; imshow(BW3_647);
%BW3_647 = imdilate(BW2_647, ones(1,1)); figure; imshow(BW3_647);

BW_488diff647 = BW3_488 .* BW3_647; %figure; imshow(BW_488diff647);
BW_488diff647_dil = imdilate(BW_488diff647, ones(5,5)); %figure; imshow(BW_488diff647_dil);

BW_488_common = BW3_488 .* BW_488diff647_dil; %figure; imshow(BW_488_common);
BW_647_common = BW3_647 .* BW_488diff647_dil; %figure; imshow(BW_647_common);
%figure; imshowpair(BW_488_common, BW_647_common)

BW_647diff488 = BW3_647 .* BW3_488; %figure; imshow(BW_647diff488);

I_488_masked = double(I_488) .* BW_488_common; %figure; imagesc(I_488_masked);
I_647_masked = double(I_647) .* BW_647_common; %figure; imagesc(I_647_masked);

p_488 = FastPeakFind(I_488_masked, 0, ones(5,5));
p_647 = FastPeakFind(I_647_masked, 0, ones(5,5));

%             figure; imagesc(I_488_masked); axis equal; hold on
%             plot(p_488(1:2:end),p_488(2:2:end),'c+')
%             plot(p_647(1:2:end),p_647(2:2:end),'m+')
% 
%             figure; imagesc(I_647_masked); axis equal; hold on
%             plot(p_488(1:2:end),p_488(2:2:end),'c+')
%             plot(p_647(1:2:end),p_647(2:2:end),'m+')

data_points488(:,1) = p_488(1:2:end-1);
data_points488(:,2) = p_488(2:2:end);

data_points488_sorted = sortrows(data_points488, [1 2]);

data_points647(:,1) = p_647(1:2:end);
data_points647(:,2) = p_647(2:2:end);

data_points647_sorted = sortrows(data_points647, [1 2]);

l488 = length(data_points488_sorted);
l647 = length(data_points647_sorted);

diff_len = abs(l488 - l647);

if l488 > l647
    data_points488_sorted = data_points488_sorted(1:end-diff_len,:);  
elseif l488 < l647
    data_points647_sorted = data_points647_sorted(1:end-diff_len,:);
end

data_diff = data_points488_sorted - data_points647_sorted;
err = abs(data_diff) < 3;

index = 1;
index_u = 1;
for kk = 1:length(data_diff)
    if err(kk,2) == 1
        data_points488_cor(index,:) = data_points488_sorted(kk,:);
        data_points647_cor(index,:) = data_points647_sorted(kk,:);
        index = index + 1;
        
    else
        points_unaligned488(index_u,:) = data_points488_sorted(kk,:);
        points_unaligned647(index_u,:) = data_points647_sorted(kk,:);
        index_u = index_u + 1;
    end
end

points_aligned488 = sortrows(points_unaligned488, 2); 
points_aligned647 = sortrows(points_unaligned647, 2);

data_diff_alig = points_aligned488 - points_aligned647;
err_2 = abs(data_diff_alig) < 2;

index_fix = 1;
for kk = 1:length(data_diff_alig)
    if err_2(kk,2) == 1
        points_aligned488_cor(index_fix,:) = points_aligned488(kk,:);
        points_aligned647_cor(index_fix,:) = points_aligned647(kk,:);
        index_fix = index_fix + 1;
    end
end

data_points488_cor = [data_points488_cor; points_aligned488_cor];
data_points647_cor = [data_points647_cor; points_aligned647_cor];

data_diff_cor = data_points488_cor - data_points647_cor;

n = length(data_points488_cor); % number of grid holes to be used for registration (x,y coordinates)
tform = fitgeotrans(data_points647_cor, data_points488_cor, 'lwm', n);

end
%
% figure; imshowpair(Grid_488_raw*10, Grid_647_reg*4, 'ColorChannels', 'green-magenta'); hold on
% plot(sorted_488(:,1),sorted_488(:,2),'cx', 'LineWidth',1, 'MarkerSize',15); hold on;
% plot(sorted_647(:,1),sorted_647(:,2),'mx', 'LineWidth',1, 'MarkerSize',15); hold off;
%
% Grid_647_reg_reg = imwarp(Grid_647_reg, tform, 'OutputView', imref2d(I));
% points647reg_raw = detectHarrisFeatures(Grid_647_reg_reg, 'FilterSize', f);
%     points647reg = sortrows(points647reg_raw.Location, [-2 1]);
% sorted_647reg = zeros(length(points647reg), 2);
%
%
%     for ff = 0:(length(points488)/gg)-1
%         sorted_647reg(1+(gg*ff):gg+(gg*ff),:) = sortrows(points647reg(1+(gg*ff):gg+(gg*ff), :), -1);
%     end
%
% figure; imshowpair(Grid_488_raw*10, Grid_647_reg_reg*4, 'ColorChannels', 'green-magenta'); hold on
% plot(sorted_488(:,1),sorted_488(:,2),'cx', 'LineWidth',1, 'MarkerSize',15); hold on;
% plot(sorted_647reg(:,1),sorted_647reg(:,2),'mx', 'LineWidth',1, 'MarkerSize',15); hold off;

%%% Notes
% [mp,fp] = cpselect(T647_Cor{1,1}{1,1}*100,I_488*50,'Wait',true);
% tform = fitgeotrans(mp, fp, 'lwm', 35);
% AA = imwarp(T647_Cor{1,1}{1,1}, tform, 'OutputView', imref2d(I));

    
    
