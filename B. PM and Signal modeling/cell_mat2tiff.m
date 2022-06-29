function cell_mat2tiff(path2save, data)
file_name = inputname(2);

    if ~iscell(data)
        TIFF = uint16(data);
        
    else
        [size1, size2] = size(data{1,1});
        lenght_data = length(data);
        TIFF = zeros(size1, size2, lenght_data);
        
        if strfind(file_name, '488') >= 1 
             TIFF = uint16(TIFF);
        elseif strfind(file_name, '647') >= 1
            TIFF = uint16(TIFF);
        end
        
        %data converter to be bfsave friendlly        
        for ii = 1:lenght_data
            TIFF(:,:,ii) = data{ii,1};
        end
    end
bfsave(TIFF, path2save);
end