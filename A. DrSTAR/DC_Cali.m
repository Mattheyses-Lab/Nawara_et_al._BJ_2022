%%% Function created by Nawara T. (Mattheyses lab - 05/20/2022) compatible
%%% with MATLAB R2020b - %This is a MATLAB version of Cairn Image Splitter
%%%(but only the calibration file generator) it is a hybrid of code and app programing

%This work is licensed under the Creative Commons Attribution 4.0
%International License. To view a copy of this license, visit
%http://creativecommons.org/licenses/by/4.0/ or send a letter to Creative
%Commons, PO Box 1866, Mountain View, CA 94042, USA.

function [Calibration_File, Grid_488_raw, Grid_647_reg] = DC_Cali(yourpath)

ContentInFold = dir(yourpath);
for i = 1:length(ContentInFold)
    if strfind(ContentInFold(i).name, 'DC_488_Cali') >= 1
        raw_grid_488 = bfopen(fullfile(yourpath,ContentInFold(i).name));
        avg_grid_488 = raw_grid_488{1, 1}(:,1);
        
        avg_grid_raw_488 = cat(3,avg_grid_488{:});
        avg_grid_488 = mean(avg_grid_raw_488,3);
        
    elseif strfind(ContentInFold(i).name, 'DC_647_Cali') >= 1
        raw_grid_647 = bfopen(fullfile(yourpath,ContentInFold(i).name));
        avg_grid_647 = raw_grid_647{1, 1}(:,1);
        
        avg_grid_raw_647 = cat(3,avg_grid_647{:});
        avg_grid_647 = mean(avg_grid_raw_647,3);
    end
end


%Pick 1st ROI then move the same ROI into the matching GRID and then
%press OK
draw_ROI = 1;
while draw_ROI == 1
    figure(1)
    imagesc(avg_grid_488, [min(avg_grid_488, [],'all'), max(avg_grid_488, [],'all')*0.3])
    axis equal
    roi1 = drawrectangle;
    ROI_488 = roi1.Position;
    imagesc(avg_grid_647, [min(avg_grid_647, [],'all'), max(avg_grid_647, [],'all')*0.3])
    axis equal
    roi2 = drawrectangle('Position',ROI_488,'StripeColor','r','InteractionsAllowed', 'translate');
    choice = menu('Press OK after moving the red ROI','OK','Redraw');
    ROI_647 = roi2.Position;
    if choice==2 || choice==0
        fprintf('Redraw the ROI');
    else
        draw_ROI = 0;
        ROI_488 = round(ROI_488);
        ROI_647 = round(ROI_647);
        Grid_488_raw = avg_grid_488(ROI_488(2):(ROI_488(2)+ROI_488(4)),ROI_488(1):(ROI_488(1)+ROI_488(3)));
        Grid_647_raw = avg_grid_647(ROI_647(2):(ROI_647(2)+ROI_647(4)),ROI_647(1):(ROI_647(1)+ROI_647(3)));
        save([yourpath, '/', 'Corrections', '/','Grid_488_raw.mat'],'Grid_488_raw');
        save([yourpath, '/', 'Corrections', '/','Grid_647_raw.mat'],'Grid_647_raw');
        assignin('base','Grid_488_raw', Grid_488_raw)
        assignin('base','Grid_647_raw', Grid_647_raw)
        % Align the grid, it supports declimals but has to be typed
        Cali_app
    end
end

choice = menu('Press OK after aligment','OK');
X = evalin('base','X');
Y = evalin('base','Y');
if isequal(X, [])
    X = 0;
end

if isequal(Y, [])
    Y = 0;
end

ROI_647(1) = ROI_647(1) - X;
ROI_647(2) = ROI_647(2) - Y;

Grid_647_raw = avg_grid_647(ROI_647(2):(ROI_647(2)+ROI_647(4)),ROI_647(1):(ROI_647(1)+ROI_647(3)));
figure(2)
imshowpair(Grid_488_raw,Grid_647_raw);
Grid_647_reg = Grid_647_raw;


delete([yourpath, '/', 'Corrections', '/','Grid_647_raw.mat'])
delete([yourpath, '/', 'Corrections', '/','Grid_488_raw.mat'])
%save([yourpath, '/', 'Corrections', '/','Grid_647_reg.mat'], 'Grid_647_reg');
Calibration_File = [ROI_488; ROI_647];
%calibation file col1 is for 488 = X, Y, then the amount you need to add to X
%tohave botm of ROI and then same but for Y, col2 is for 647
path488 = [yourpath, '/', 'Corrections', '/','Grid_488_raw.tif'];
path647 = [yourpath, '/', 'Corrections', '/','Grid_647_reg.tif'];
cell_mat2tiff(path488, Grid_488_raw)
cell_mat2tiff(path647, Grid_647_reg)
save([yourpath, '/', 'Corrections', '/','Calibration_File.mat'],'Calibration_File');
close all
end


