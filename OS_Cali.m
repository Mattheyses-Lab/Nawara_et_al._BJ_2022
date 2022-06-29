%%% Function created by Nawara T. (Mattheyses lab - 05/20/2022) compatible
%%% with MATLAB R2020b - %This is a MATLAB version of Cairn Image Splitter 
%%%(but only the calibration file generator) it is a hybrid of code and app programing

%This work is licensed under the Creative Commons Attribution 4.0
%International License. To view a copy of this license, visit
%http://creativecommons.org/licenses/by/4.0/ or send a letter to Creative
%Commons, PO Box 1866, Mountain View, CA 94042, USA.

function [OS_Calibration_File, Grid_488_raw, Grid_647_reg] = OS_Cali(yourpath)

ContentInFold = dir(yourpath);
    for i = 1:length(ContentInFold)
        if strfind(ContentInFold(i).name, 'OS_Cali') >= 1  
            raw_grid = bfopen(fullfile(yourpath,ContentInFold(i).name));
            avg_grid = raw_grid{1, 1}(:,1);
            
            avg_grid_raw = cat(3,avg_grid{:});
            avg_grid = mean(avg_grid_raw,3);
            
            %Pick 1st ROI then move the same ROI into the matching GRID and then
            %press OK
            figure(1)
            imagesc(avg_grid)
            axis equal
            roi1 = drawrectangle;
            roi2 = drawrectangle('Position',[roi1.Position],'StripeColor','r','InteractionsAllowed', 'translate');
            choice = menu('Press OK after moving the red ROI','OK','Cancel');
                if choice==2 || choice==0
                   return;
                else
                    roi1.Position = round(roi1.Position);
                    roi2.Position = round(roi2.Position);
                    Grid_488_raw = avg_grid(roi1.Position(2):(roi1.Position(2)+roi1.Position(4)),roi1.Position(1):(roi1.Position(1)+roi1.Position(3)));
                    Grid_647_raw = avg_grid(roi2.Position(2):(roi2.Position(2)+roi2.Position(4)),roi2.Position(1):(roi2.Position(1)+roi2.Position(3)));
                    save([yourpath, '/', 'Corrections', '/','Grid_488_raw.mat'],'Grid_488_raw');
                    save([yourpath, '/', 'Corrections', '/','Grid_647_raw.mat'],'Grid_647_raw');
                    assignin('base','Grid_488_raw', Grid_488_raw)
                    assignin('base','Grid_647_raw', Grid_647_raw)
                    % Align the grid, it supports declimals but has to be typed   
                    OS_Cali_app        
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

             roi2.Position(1) = roi2.Position(1) - X;
             roi2.Position(2) = roi2.Position(2) - Y;

             Grid_647_raw = avg_grid(roi2.Position(2):(roi2.Position(2)+roi2.Position(4)),roi2.Position(1):(roi2.Position(1)+roi2.Position(3)));
             figure(2)
             imshowpair(Grid_488_raw,Grid_647_raw);
             Grid_647_reg = Grid_647_raw;


             delete([yourpath, '/', 'Corrections', '/','Grid_647_raw.mat'])
             delete([yourpath, '/', 'Corrections', '/','Grid_488_raw.mat'])
             %save([yourpath, '/', 'Corrections', '/','Grid_647_reg.mat'], 'Grid_647_reg');
             OS_Calibration_File = [roi1.Position; roi2.Position];
             %calibation file col1 is for 488 = X, Y, then the amount you need to add to X
             %tohave botm of ROI and then same but for Y, col2 is for 647
             path488 = [yourpath, '/', 'Corrections', '/','Grid_488_raw.tif'];
             path647 = [yourpath, '/', 'Corrections', '/','Grid_647_reg.tif'];
             cell_mat2tiff(path488, Grid_488_raw)
             cell_mat2tiff(path647, Grid_647_reg)
             save([yourpath, '/', 'Corrections', '/','OS_Calibration_File.mat'],'OS_Calibration_File');
             close all
        end
    end
end

