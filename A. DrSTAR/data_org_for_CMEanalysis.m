%%% Function created by Nawara T. (Mattheyses lab - 05/20/2022) compatible
%%% with MATLAB R2020b / it organizes data post DrSTAR run to so it can be
%%% transfered directlly into CMEanalysis, if multiple condtions are
%%% included user has to group them by had as examination of inteisty and
%%% dz channel should be done prior trnasferring for CMEanalysis

%%% CRITICAL: Adjust frame_rate based on experiamntal data!!!!

%This work is licensed under the Creative Commons Attribution 4.0
%International License. To view a copy of this license, visit
%http://creativecommons.org/licenses/by/4.0/ or send a letter to Creative
%Commons, PO Box 1866, Mountain View, CA 94042, USA.

function data_org_for_CMEanalysis(yourpath)

dd = 1;
frame_rate = 0.3;

ContentInFold = dir([yourpath, '/Data_analysis']);
mkdir([yourpath, '/', 'For_CMEanalysis']);

    
    for i = 1:length(ContentInFold)
        if strfind(ContentInFold(i).name, 'STAR') >= 1
            mkdir([yourpath, '/', 'For_CMEanalysis', '/', ['Cell', num2str(dd), '_', num2str(frame_rate), 's']])
            dir_Ch1 = [yourpath, '/', 'For_CMEanalysis', '/', ['Cell', num2str(dd), '_', num2str(frame_rate), 's'], '/', 'Ch1'];
            dir_Ch2 = [yourpath, '/', 'For_CMEanalysis', '/', ['Cell', num2str(dd), '_', num2str(frame_rate), 's'], '/', 'Ch2'];
            dir_Ch3 = [yourpath, '/', 'For_CMEanalysis', '/', ['Cell', num2str(dd), '_', num2str(frame_rate), 's'], '/', 'Ch3'];
            mkdir(dir_Ch1)
            mkdir(dir_Ch2)
            mkdir(dir_Ch3)
            T488_Cor_Blur_path = [yourpath, '/', 'Data_analysis', '/', ContentInFold(i).name, '/', ContentInFold(i).name, '_488_Cor_Blur.tif'];
            T647_Cor_Reg_path = [yourpath, '/', 'Data_analysis', '/', ContentInFold(i).name, '/', ContentInFold(i).name, '_647_Cor_Reg.tif'];
            dz_channel_path = [yourpath, '/', 'Data_analysis', '/', ContentInFold(i).name, '/', ContentInFold(i).name, '_dz_channel.tif'];
            movefile(T488_Cor_Blur_path, dir_Ch1);
            movefile(T647_Cor_Reg_path, dir_Ch2);
            movefile(dz_channel_path, dir_Ch3);
            dd = dd + 1;
        end
    end
    writematrix([], [yourpath, '/', 'Corrections', '/','organization_done.txt']); %added 032122 for automation if somehting gees wrong
end