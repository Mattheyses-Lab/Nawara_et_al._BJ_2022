% % % Created by Nawara T. Matttheyses Lab @UAB 03.02.20
% % % IMREG is a fucntion that allows to register multichannel image stack,
% % % based on the calibration grid. it needs 7XXX types of imputs:
% % % 1. Path to data folder
% % % 2. Referance grid picture [Grid_488_raw];
% % % 2. Aligned grid picture [Grid_647_reg];
% % % Program requires instalation of bio-formats tool box details:
% % % https://docs.openmicroscopy.org/bio-formats/6.4.0/users/matlab/index.html
% % % Requires Computer Vision Toolbox
% % % For any qestions email tomasz.nawara.95@gmail.com

%This work is licensed under the Creative Commons Attribution 4.0
%International License. To view a copy of this license, visit
%http://creativecommons.org/licenses/by/4.0/ or send a letter to Creative
%Commons, PO Box 1866, Mountain View, CA 94042, USA.

%[last edited:05.19.22]
%  Fixed error in how detected points were saved in matrix [04.26.20]
%  Fixed image display [04.27.20]
%  Added the registration conformation window (figure 1) [05.05.20]
%  Added Registration of signle frame images [07.04.20]
%  Added batch processing [10.09.20]
%  Added batch processing for single images [11.10.20]
%  Display for registration confomration [04.11.22]
%  Changed registration for loop to 3D matrix [05.19.22]
%  Added data based registtration and manaul registration [06.29.22]

function IMGREG_V2(yourpath, Grid_488_raw, Grid_647_reg)

stop_it_please = 0; % while loop for registration

ContentInFold = dir(yourpath);

% registration mask generation for every cell
for i = 1:length(ContentInFold)
    if strfind(ContentInFold(i).name, 'STAR') >= 1
        T488_cor_path = [yourpath, '/', 'Data_analysis', '/', ContentInFold(i).name(1:end-4), '/', ContentInFold(i).name(1:end-4), '_488_Cor.tif'];
        T647_cor_path = [yourpath, '/', 'Data_analysis', '/', ContentInFold(i).name(1:end-4), '/', ContentInFold(i).name(1:end-4), '_647_Cor.tif'];
        fn647 = [yourpath, '/', 'Data_analysis', '/', ContentInFold(i).name(1:end-4), '/', ContentInFold(i).name(1:end-4), '_647_Cor_Reg.tif'];
        tform_path = [yourpath, '/', 'Data_analysis', '/', ContentInFold(i).name(1:end-4), '/' ContentInFold(i).name(1:end-4), '_tform.mat'];
        fig_path = [yourpath, '/', 'Data_analysis', '/', ContentInFold(i).name(1:end-4), '/' ContentInFold(i).name(1:end-4), '_reg_fig.png'];
        
        %loads existing tform
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
            
            %trys to register based on data if fails moves on
            try
                [tform_data, ~, ~] = tform_generation_data(I_488, I_647);
                I_647_data = imwarp(I_647, tform_data, 'OutputView', imref2d(I));
            catch
                fprintf('Data registration not possible to few points')
                I_647_data = ones(size(I_647));
            end
            
            
            [tform_grid, ~, ~] = tform_generation_grid(Grid_488_raw, Grid_647_reg);
            I_647_grid = imwarp(I_647, tform_grid, 'OutputView', imref2d(I));
            
            while stop_it_please == 0
                
                figure('units','normalized','outerposition',[0 0 1 1])
                tiledlayout(2,3);
                nexttile; imshowpair(imadjust(I_488), imadjust(I_647_raw), 'ColorChannels', 'green-magenta'); title('Unregistered overlay');
                nexttile; imshowpair(imadjust(imgaussfilt(I_488,1)), imadjust(I_647_data), 'ColorChannels', 'green-magenta'); title('Data based registration');
                nexttile; imshowpair(imadjust(imgaussfilt(I_488,1)), imadjust(I_647_grid), 'ColorChannels', 'green-magenta'); title('Grid based registration');
                nexttile; imshow(double(I_647_raw)./double(I_488), [0 3]); title('Raw 647/488 ratio'); hold on; set(gca,'xticklabel',{[]}); set(gca,'yticklabel',{[]})
                nexttile; imshow(double(I_647_data)./double(imgaussfilt(I_488,1)), [0 3]); title('Data 647/488 ratio'); hold on; set(gca,'xticklabel',{[]}); set(gca,'yticklabel',{[]})
                nexttile; imshow(double(I_647_grid)./double(imgaussfilt(I_488,1)), [0 3]); title('Grid 647/488 ratio'); hold on; set(gca,'xticklabel',{[]}); set(gca,'yticklabel',{[]})
                
                %Chose the registration type D (Data based registration) or G (Grid based registration) if results are not
                %satisfactroy press M (manual XY selction)
                m = input('Chose Data/Grid/Manual registration (D/G/M):','s');
                
                if m == 'D'
                    
                    FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
                    save_tfrom(tform_path, tform_data, fig_path, FigList)
                    close all
                    stop_it_please = 1;
                    
                elseif m == 'G'
                    
                    FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
                    save_tfrom(tform_path, tform_grid, fig_path, FigList)
                    close all
                    stop_it_please = 1;
                    
                elseif m == 'M'
                    %chose maching XYs one by one
                    close all
                    [mp,fp] = cpselect(imadjust(I_647_raw),imadjust(I_488),'Wait',true);
                    n = length(mp);
                    tform_manual = fitgeotrans(mp, fp, 'lwm', n);
                    I647_manual = imwarp(I_647_raw, tform_manual, 'OutputView', imref2d(I));
                    
                    tiledlayout(2,1);
                    nexttile; imshowpair(imadjust(imgaussfilt(I_488,1)), imadjust(I647_manual), 'ColorChannels', 'green-magenta'); title('Manualy registered');
                    nexttile; imshow(double(I647_manual)./double(imgaussfilt(I_488,1)), [0 3]); title('Manual 647/488 ratio'); hold on; set(gca,'xticklabel',{[]}); set(gca,'yticklabel',{[]})
                    
                    m = input('Do you want  to continue (Y/N):','s');
                    
                    if m == 'N'
                        stop_it_please = 0;
                        close all
                    elseif m == 'Y'
                        FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
                        save_tfrom(tform_path, tform_manual, fig_path, FigList)
                        stop_it_please = 1;
                    end
                end
            end
            stop_it_please = 0;
        end
    end
end

%image regsitraiton
parfor i = 1:length(ContentInFold)
    if strfind(ContentInFold(i).name, 'STAR') >= 1
        T647_cor_path = [yourpath, '/', 'Data_analysis', '/', ContentInFold(i).name(1:end-4), '/', ContentInFold(i).name(1:end-4), '_647_Cor.tif'];
        fn647 = [yourpath, '/', 'Data_analysis', '/', ContentInFold(i).name(1:end-4), '/', ContentInFold(i).name(1:end-4), '_647_Cor_Reg.tif']; %%%%
        tform_path = [yourpath, '/', 'Data_analysis', '/', ContentInFold(i).name(1:end-4), '/' ContentInFold(i).name(1:end-4), '_tform.mat'];
        tform = load(tform_path);
        tform = tform.tform;
        
        if ~exist(fn647, 'file') %%%%
            
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
        end %%$%%
    end
end
writematrix([], [yourpath, '/', 'Corrections', '/','647_Registration_done.txt']); %added 032122 for automation if somehting gees wrong
fprintf('<3 Registration completed <3 \n');
end

function save_tfrom(tform_path, tform, fig_path, FigList)
save(tform_path, 'tform');
try
    saveas(FigList, fig_path)
catch
    fprintf('Registration figure was closed and will not be saved')
end
end

%generation of data based tform
function [tform, p_488, p_647] = tform_generation_data(I_488, I_647)

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
BW2_488_size = bwpropfilt(BW_488, 'Area', [30 100]); %figure; imshow(BW2_488_size);
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
BW2_647_size = bwpropfilt(BW_647, 'Area', [30 100]); %figure; imshow(BW2_647_size);
BW3_647 = bwpropfilt(BW2_647_size, "Eccentricity", [0 0.55]); %figure; imshow(BW3_647);
%BW3_647 = imdilate(BW2_647, ones(1,1)); figure; imshow(BW3_647);

BW_488diff647 = BW3_488 .* BW3_647; %figure; imshow(BW_488diff647);
BW_488diff647_dil = imdilate(BW_488diff647, ones(5,5)); %figure; imshow(BW_488diff647_dil);

BW_488_common = BW3_488 .* BW_488diff647_dil; %figure; imshow(BW_488_common);
BW_647_common = BW3_647 .* BW_488diff647_dil; %figure; imshow(BW_647_common);
% figure; imshowpair(BW_488_common, BW_647_common)

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


data_diff_cor = data_points488_cor - data_points647_cor;


n = length(data_points488_cor); % number of grid holes to be used for registration (x,y coordinates)
tform = fitgeotrans(data_points647_cor, data_points488_cor, 'lwm', n);

end

%generation of grid based tform
function [tform, sorted_488, sorted_647] = tform_generation_grid(Grid_488_raw, Grid_647_reg)

sigma = 3;
f = 11;

%if use_fixed_detection == 0
%Apply Gausian smotching to pictures
I488 = imgaussfilt(Grid_488_raw, sigma);
I647 = imgaussfilt(Grid_647_reg, sigma);

%Use bfopen to load the data
I = size(Grid_488_raw); %size of the referance picture

%Detect features
points488_raw = detectHarrisFeatures(I488, 'FilterSize', f);
points488 = sortrows(points488_raw.Location, [-2 1]);

points647_raw = detectHarrisFeatures(I647, 'FilterSize', f);
points647 = sortrows(points647_raw.Location, [-2 1]);

%     figure(1)
%     tiledlayout(2,2)
%     nexttile; imagesc(Grid_488_raw); title('Raw 488')
%     nexttile; imshow(imadjust(I488)); hold on; plot(points488_raw); title('Detection 488'); hold off;
%     nexttile; imagesc(Grid_647_reg); title('Raw 647')
%     nexttile; imshow(imadjust(I488)); hold on; plot(points488_raw); title('Detection 647'); hold off;

sorted_488 = zeros(length(points488), 2);
sorted_647 = zeros(length(points647), 2);

%Accept detection
%     m=input('Do you want to continue, Y/N:','s');
%     if m=='N'
%         close all
%         return
%
%     elseif m=='Y'
%Create registration mask
%but correct for terrible matching of
%detection grid
%if you delet numbers bigger from around xzero then you automatically
%elimnaets unaligned points :) done
for gg = 1:length(points488)
    if points488(gg,2) - points488(gg+1,2) > 8
        break
    end
end

for ff = 0:(length(points488)/gg)-1
    sorted_488(1+(gg*ff):gg+(gg*ff),:) = sortrows(points488(1+(gg*ff):gg+(gg*ff), :), -1);
    sorted_647(1+(gg*ff):gg+(gg*ff),:) = sortrows(points647(1+(gg*ff):gg+(gg*ff), :), -1);
end

n = length(sorted_488); % number of grid holes to be used for registration (x,y coordinates)
tform = fitgeotrans(sorted_647, sorted_488, 'lwm', n);
%end
end


%Notes
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