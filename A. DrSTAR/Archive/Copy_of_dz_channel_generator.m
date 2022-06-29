function Copy_of_dz_channel_generator(yourpath)

gamma = 0.0024; %gamma value calulted form STAR calculator
%sigma = 1.01; % a blure that will event out the registration on 647 on a 488 channel %%%ADD Automation here for sigma calcultion


    % Get all the subfolders
    ContentInFold = dir(yourpath);
    
    parfor i = 1:length(ContentInFold)
        if strfind(ContentInFold(i).name, 'Cell') >= 1
            T488_cor_path = [yourpath, '/', ContentInFold(i).name, '/', 'Ch1'];
            ContentInSubFold = dir(T488_cor_path);
            T488 = bfopen([ContentInSubFold(3).folder, '/', ContentInSubFold(3).name]);
            
            T647_cor_path = [yourpath, '/', ContentInFold(i).name, '/', 'Ch2'];
            ContentInSubFold = dir(T647_cor_path);
            T647 = bfopen([ContentInSubFold(3).folder, '/', ContentInSubFold(3).name]);
                          
           
            dz_channel = cell(size(T488{1,1},1),1);
            for ii = 1:length(dz_channel)
                dz_channel{ii,1} = zeros(size(T488{1,1}{1,1}));
            end
            
            ratio_488 = dz_channel;
            ratio_488_10 = dz_channel{1,1};
            ratio_647 = dz_channel;
            ratio_647_10 = dz_channel{1,1};
            
            
            for z = 1:10
                ratio_488_10 = ratio_488_10 + double(T488{1,1}{z,1});
                ratio_647_10 = ratio_647_10 + double(T647{1,1}{z,1});
            end
            
            ratio_488_10 = ratio_488_10/10;
            ratio_647_10 = ratio_647_10/10;
            
            for k = 1:length(T488{1,1})   
                ratio_488{k,1} = double(T488{1,1}{k, 1}) ./ ratio_488_10;
                ratio_647{k,1} = double(T647{1,1}{k, 1}) ./ ratio_647_10;
                ratio_647{k,1} = [ratio_647{k,1}(2:end, 1:end); zeros(1, size(ratio_647{k,1},2))];
                dz_channel{k,1} =  log(ratio_647{k,1} ./ ratio_488{k,1})/gamma;
            end
            
            
            fndz = [yourpath, '/', ContentInFold(i).name, '/', 'Ch3', '/', 'new_dz_channel.tif'];
            cell_mat2tiff(fndz, dz_channel);
        end
    end
end