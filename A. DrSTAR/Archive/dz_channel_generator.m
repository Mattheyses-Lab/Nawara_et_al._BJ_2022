function dz_channel_generator(yourpath, y)

gamma = gamma_calcultor(y, yourpath); %gamma value calulted form STAR calculator

try
    sigma = evalin('base','Sigma_cor');
catch
    sigma = 1.01; % a blure that will event out the registration on 647 on a 488 channel %%%ADD Automation here for sigma calcultion
end

% Get all the subfolders
ContentInFold = dir(yourpath);

for i = 1:length(ContentInFold)
    if strfind(ContentInFold(i).name, 'STAR') >= 1
        T488_cor_path = [yourpath, '/', 'Data_analysis', '/', ContentInFold(i).name(1:end-4), '/', ContentInFold(i).name(1:end-4), '_488_Cor.tif'];
        T647_cor_reg_path = [yourpath, '/', 'Data_analysis', '/', ContentInFold(i).name(1:end-4), '/', ContentInFold(i).name(1:end-4), '_647_Cor_Reg.tif'];
        T488_Cor = bfopen(T488_cor_path);
        
        T488_Cor_Blur = cell(size(T488_Cor{1,1},1),1);
        for ii = 1:size(T488_Cor_Blur{1,1},1)
            T488_Cor_Blur{ii,1} = zeros(size(T488_Cor{1,1}{1,1}));
        end
        
        for k = 1:size(T488_Cor{1,1},1) %Generate blured 488 to match 647 intrploation
            T488_Cor_Blur{k, 1} = imgaussfilt(T488_Cor{1,1}{k, 1}, sigma);
        end

        T488_Cor = [];
        
        dz_channel = cell(size(T488_Cor_Blur));
        for ii = 1:size(T488_Cor_Blur,1)
            dz_channel{ii,1} = zeros(size(T488_Cor_Blur{1,1}));
        end
        
        ratio_488 = dz_channel;
        %ratio_488_10 = dz_channel{1,1};
        ratio_647 = dz_channel;
        %ratio_647_10 = dz_channel{1,1};
        
        T647_Cor_Reg2open = bfopen(T647_cor_reg_path);
        
        T647_Cor_Reg = cell(size(T647_Cor_Reg2open{1,1},1),1);
        for ii = 1:size(T647_Cor_Reg,1)
            T647_Cor_Reg{ii,1} = T647_Cor_Reg2open{1,1}{ii,1};
        end
        
        [PM_only_488, Cell_mask] = Signal_removal(T488_Cor_Blur, 488, []);
        PMfn488 = [yourpath, '/', 'Data_analysis', '/', ContentInFold(i).name(1:end-4), '/', ContentInFold(i).name(1:end-4), '_488_PM.tif'];
        cell_mat2tiff(PMfn488, PM_only_488);

        [PM_only_647, ~] = Signal_removal(T647_Cor_Reg, 647, Cell_mask);
        PMfn647 = [yourpath, '/', 'Data_analysis', '/', ContentInFold(i).name(1:end-4), '/', ContentInFold(i).name(1:end-4), '_647_PM.tif'];
        cell_mat2tiff(PMfn647, PM_only_647);
        
        path_cell_mask = [yourpath, '/', 'Data_analysis', '/', ContentInFold(i).name(1:end-4), '/', ContentInFold(i).name(1:end-4), '_cell_mask.tif'];
        cell_mat2tiff(path_cell_mask, Cell_mask);

        
        %%%%%%%%% ADD NEW DZ channel generator
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
        
        for k = 1:size(T488_Cor_Blur,1)
            ratio_488{k,1} = double(T488_Cor_Blur{k, 1}) ./ double(PM_only_488{k, 1});
            ratio_647{k,1} = double(T647_Cor_Reg{k, 1}) ./ double(PM_only_647{k, 1});
            dz_channel{k,1} =  log(ratio_647{k,1} ./ ratio_488{k,1})/gamma;
        end
        
        fn488 = [yourpath, '/', 'Data_analysis', '/', ContentInFold(i).name(1:end-4), '/', ContentInFold(i).name(1:end-4), '_488_Cor_Blur.tif'];
        cell_mat2tiff(fn488, T488_Cor_Blur);
        
        T488_Cor_Blur = [];
        T647_Cor_Reg = [];
        PM_only_488 = [];  
        PM_only_647 = []; 
        
        fndz = [yourpath, '/', 'Data_analysis', '/', ContentInFold(i).name(1:end-4), '/', ContentInFold(i).name(1:end-4), '_dz_channel.tif'];
        cell_mat2tiff(fndz, dz_channel);
    end
end
 writematrix([], [yourpath, '/', 'Corrections', '/','dz_generation_done.txt']); %added 032122 for automation if somehting gees wrong
end