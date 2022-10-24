%%% Function created by Nawara T. (Mattheyses lab - 05/20/2022) compatible
%%% with MATLAB R2020b / Function Runs blur correction on 488, creates
%%% background plasma mabrane pictrue for both channels, creates ratio and
%%% ultimatelly the dZ channel adjust sigma accordingly with your averaged
%%% interpolation value

%This work is licensed under the Creative Commons Attribution 4.0
%International License. To view a copy of this license, visit
%http://creativecommons.org/licenses/by/4.0/ or send a letter to Creative
%Commons, PO Box 1866, Mountain View, CA 94042, USA.

function dz_channel_generator(yourpath, y)

%gamma value calulted form STAR calculator
gamma = gamma_calculator(y);

% a blur that will event out the registration on 647 on a 488 channel
try
    sigma = evalin('base','Sigma_cor');
catch
    sigma = 1.01; 
end

% Get all the subfolders
ContentInFold = dir([yourpath, '/Data_analysis']);

parfor i = 1:length(ContentInFold)
    if strfind(ContentInFold(i).name, 'STAR') >= 1
        T488_cor_path = [yourpath, '/', 'Data_analysis', '/', ContentInFold(i).name, '/', ContentInFold(i).name, '_488_Cor.tif'];
        T647_cor_reg_path = [yourpath, '/', 'Data_analysis', '/', ContentInFold(i).name, '/', ContentInFold(i).name, '_647_Cor_Reg.tif'];
        fn488 = [yourpath, '/', 'Data_analysis', '/', ContentInFold(i).name, '/', ContentInFold(i).name, '_488_Cor_Blur.tif'];
        PMfn488 = [yourpath, '/', 'Data_analysis', '/', ContentInFold(i).name, '/', ContentInFold(i).name, '_488_PM.tif'];
        PMfn488_signal = [yourpath, '/', 'Data_analysis', '/', ContentInFold(i).name, '/', ContentInFold(i).name, '_488_Signal_mask.tif'];
        PMfn647 = [yourpath, '/', 'Data_analysis', '/', ContentInFold(i).name, '/', ContentInFold(i).name, '_647_PM.tif'];
        fndz = [yourpath, '/', 'Data_analysis', '/', ContentInFold(i).name, '/', ContentInFold(i).name, '_dz_channel.tif'];
        path_cell_mask = [yourpath, '/', 'Data_analysis', '/', ContentInFold(i).name, '/', ContentInFold(i).name, '_cell_mask.tif'];
        
        if ~exist(fndz, 'file')
            if exist(fn488, 'file')
                T488_Cor_Blur = bfopen(fn488);
                T488_Cor_Blur = T488_Cor_Blur{1,1}(:,1);
            else
                T488_Cor = bfopen(T488_cor_path);
                T488_Cor_Blur = cell(size(T488_Cor{1,1},1),1);
                for ii = 1:size(T488_Cor_Blur{1,1},1)
                    T488_Cor_Blur{ii,1} = zeros(size(T488_Cor{1,1}{1,1}));
                end
                
                %Generate blured 488 to match 647 intrploation
                for k = 1:size(T488_Cor{1,1},1) 
                    T488_Cor_Blur{k, 1} = imgaussfilt(T488_Cor{1,1}{k, 1}, sigma);
                end
                T488_Cor = [];
            end
            
            dz_channel = cell(size(T488_Cor_Blur));
            for ii = 1:size(T488_Cor_Blur,1)
                dz_channel{ii,1} = zeros(size(T488_Cor_Blur{1,1}));
            end
            
            ratio_488 = dz_channel;
            %ratio_488_10 = dz_channel{1,1};
            ratio_647 = dz_channel;
            %ratio_647_10 = dz_channel{1,1};
            
            T647_Cor_Reg2open = bfopen(T647_cor_reg_path);
            T647_Cor_Reg = T647_Cor_Reg2open{1,1}(:,1);
            
            %         T647_Cor_Reg = cell(size(T647_Cor_Reg2open{1,1},1),1);
            %         for ii = 1:size(T647_Cor_Reg,1)
            %             T647_Cor_Reg{ii,1} = T647_Cor_Reg2open{1,1}{ii,1};
            %         end
            
            %Generation of plasma membrane picture for 488 emission
            if exist(PMfn488, 'file')
                PM_only_488 = bfopen(PMfn488);
                PM_only_488 = PM_only_488{1,1}(:,1);
                Cell_mask = bfopen(path_cell_mask);
                Cell_mask = Cell_mask{1, 1}{1, 1};  
            else
                [PM_only_488, Cell_mask, Signal_mask] = Signal_removal(T488_Cor_Blur, 488, []);
                %PMfn488 = [yourpath, '/', 'Data_analysis', '/', ContentInFold(i).name(1:end-4), '/', ContentInFold(i).name(1:end-4), '_488_PM.tif'];
                cell_mat2tiff(path_cell_mask, Cell_mask);
                cell_mat2tiff(PMfn488, PM_only_488);
                cell_mat2tiff(PMfn488_signal, Signal_mask);
            end
            
            %Generation of plasma membrane picture for 647 emission
            if exist(PMfn647, 'file')
                PM_only_647 = bfopen(PMfn647);
                PM_only_647 = PM_only_647{1,1}(:,1);
            else
                [PM_only_647, ~, ~] = Signal_removal(T647_Cor_Reg, 647, Cell_mask);
                %PMfn647 = [yourpath, '/', 'Data_analysis', '/', ContentInFold(i).name(1:end-4), '/', ContentInFold(i).name(1:end-4), '_647_PM.tif'];
                cell_mat2tiff(PMfn647, PM_only_647);
            end
                        
            %%%% DEPRACTED
            %             for z = 1:10
            %                 ratio_488_10 = ratio_488_10 + double(T488_Cor_Blur{z,1});
            %                 ratio_647_10 = ratio_647_10 + double(T647_Cor_Reg{z,1});
            %             end
            %
            %             ratio_488_10 = ratio_488_10/10;
            %             ratio_647_10 = ratio_647_10/10;
            %
            %             for k = 1:length(T488_Cor_Blur)
            %                 ratio_488{k,1} = double(T488_Cor_Blur{k, 1}) ./ ratio_488_10;
            %                 ratio_647{k,1} = double(T647_Cor_Reg{k, 1}) ./ ratio_647_10;
            %                 dz_channel{k,1} =  log(ratio_647{k,1} ./ ratio_488{k,1})/gamma;
            %             end
            %%%%
            
            %Generation of dz channel
            for k = 1:size(T488_Cor_Blur,1)
                ratio_488{k,1} = double(T488_Cor_Blur{k, 1}) ./ double(PM_only_488{k, 1});
                ratio_647{k,1} = double(T647_Cor_Reg{k, 1}) ./ double(PM_only_647{k, 1});
                dz_channel{k,1} =  log(ratio_647{k,1} ./ ratio_488{k,1})/gamma;
            end
            
            %fn488 = [yourpath, '/', 'Data_analysis', '/', ContentInFold(i).name(1:end-4), '/', ContentInFold(i).name(1:end-4), '_488_Cor_Blur.tif'];
            if ~exist(fn488, 'file')
                cell_mat2tiff(fn488, T488_Cor_Blur);
            end
            
            T488_Cor_Blur = [];
            T647_Cor_Reg = [];
            PM_only_488 = [];
            PM_only_647 = [];
            
            %fndz = [yourpath, '/', 'Data_analysis', '/', ContentInFold(i).name(1:end-4), '/', ContentInFold(i).name(1:end-4), '_dz_channel.tif'];
            cell_mat2tiff(fndz, dz_channel);
        end
    end
end
writematrix([], [yourpath, '/', 'Corrections', '/','dz_generation_done.txt']); %added 032122 for automation if somehting gees wrong
end