%%% Function created by Nawara T. (Mattheyses lab - 05/20/2022) compatible
%%% with MATLAB R2020b / this is the fuction that saves images as TIFF
%%% files either 16-bit for intensity or 32-bit for dZ channel

%This work is licensed under the Creative Commons Attribution 4.0
%International License. To view a copy of this license, visit
%http://creativecommons.org/licenses/by/4.0/ or send a letter to Creative
%Commons, PO Box 1866, Mountain View, CA 94042, USA.

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