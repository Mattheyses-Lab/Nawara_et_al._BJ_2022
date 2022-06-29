function dz_channel_generator_647_translation(yourpath, y)

gamma = gamma_calcultor(y, yourpath); %gamma value calulted form STAR calculator


    % Get all the subfolders
    ContentInFold = dir(yourpath);
    
    for i = 1:length(ContentInFold)
        if strfind(ContentInFold(i).name, 'STAR') >= 1
            T488_cor_path = [yourpath, '/', 'Data_analysis', '/', ContentInFold(i).name(1:end-4), '/', ContentInFold(i).name(1:end-4), '_488_Cor_Blur.tif'];
            T647_cor_reg_path = [yourpath, '/', 'Data_analysis', '/', ContentInFold(i).name(1:end-4), '/', ContentInFold(i).name(1:end-4), '_647_Cor_Reg.tif'];
            T488_Cor = bfopen(T488_cor_path);
            
            T488_Cor_Blur = T488_Cor{1,1};
            
            T488_Cor = [];

            dz_channel = cell(size(T488_Cor_Blur));
            for ii = 1:length(T488_Cor_Blur)
                dz_channel{ii,1} = zeros(size(T488_Cor_Blur{1,1}));
            end
            
            ratio_488 = dz_channel;
            ratio_488_10 = dz_channel{1,1};
            ratio_647 = dz_channel;
            ratio_647_10 = dz_channel{1,1};
            
            T647_Cor_Reg = bfopen(T647_cor_reg_path);
            
            for ii = 1:length(T647_Cor_Reg{1,1})
                T647_Cor_Reg{1,1}{ii,1} = imtranslate(T647_Cor_Reg{1,1}{ii,1},[0, 1]);
            end
            
            T647_Cor_Reg_Trans = T647_Cor_Reg{1,1};
            
            fn647_trans = [yourpath, '/', 'Data_analysis', '/', ContentInFold(i).name(1:end-4), '/', ContentInFold(i).name(1:end-4), '_647_Cor_Reg_Trans.tif'];
            cell_mat2tiff(fn647_trans, T647_Cor_Reg_Trans);
                       
            for z = 1:10
                ratio_488_10 = ratio_488_10 + double(T488_Cor_Blur{z,1});
                ratio_647_10 = ratio_647_10 + double(T647_Cor_Reg{1,1}{z,1});
            end
            
            ratio_488_10 = ratio_488_10/10;
            ratio_647_10 = ratio_647_10/10;
            
            for k = 1:length(T488_Cor_Blur)   
                ratio_488{k,1} = double(T488_Cor_Blur{k, 1}) ./ ratio_488_10;
                ratio_647{k,1} = double(T647_Cor_Reg{1,1}{k, 1}) ./ ratio_647_10;
                dz_channel{k,1} =  log(ratio_647{k,1} ./ ratio_488{k,1})/gamma;
            end
            
            T488_Cor_Blur = [];
            T647_Cor_Reg = [];
            
            fndz = [yourpath, '/', 'Data_analysis', '/', ContentInFold(i).name(1:end-4), '/', ContentInFold(i).name(1:end-4), '_dz_channel.tif'];
            cell_mat2tiff(fndz, dz_channel);
        end
    end
end