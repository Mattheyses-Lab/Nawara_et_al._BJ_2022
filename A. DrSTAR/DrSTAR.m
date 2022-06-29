%%% Function created by Nawara T. (Mattheyses lab - 05/20/2022) compatible
%%% with MATLAB R2020b
%%% Dream STAR
%%% Data spliting, corection, dz_channel generator, and data export
%%% package

%This work is licensed under the Creative Commons Attribution 4.0
%International License. To view a copy of this license, visit
%http://creativecommons.org/licenses/by/4.0/ or send a letter to Creative
%Commons, PO Box 1866, Mountain View, CA 94042, USA.

function DrSTAR

close all

path2code = pwd;
addpath(genpath([path2code, '\bfmatlab']));
addpath(genpath([path2code, '\image registration']));

%1 if you want the program to shutdown the computer after completion
system_shutdown = 0;
%1 if you want the data to be sorted for CMEanalysis
sort4CMEanalysis = 1;
%1 if you have list of registered XY coordiantes for IMGREG_V2

% initialization of parpool
try
    poolobj = gcp;
catch
    poolobj = parpool("threads");
end

selpath = uigetdir([], 'Select folder with the raw data');
mkdir([selpath, '/', 'Corrections']); warning('off');
mkdir([selpath, '/', 'Data_analysis']);  warning('off');

% y value is the distance form the line of the bottom of the objective
% to the laser smear on the wall (Fig S1)
if ~isfile([selpath, '/', 'Corrections', '/', 'y.mat'])
    prompt = 'What is the y value [cm]?   ';
    y = input(prompt);
    y_path = [selpath, '/', 'Corrections', '/', 'y.mat'];
    save(y_path,'y');
else
    load([selpath, '/', 'Corrections', '/', 'y.mat'])
end

%Grid splitting and Optosplit calibration
if ~isfile([selpath, '/', 'Corrections', '/', 'OS_Calibration_File.mat'])
    [OS_Calibration_File, Grid_488_raw, Grid_647_reg] = OS_Cali(selpath);
else
    r_488 = bfGetReader([selpath, '/', 'Corrections', '/', 'Grid_488_raw.tif']);
    r_647 = bfGetReader([selpath, '/', 'Corrections', '/', 'Grid_647_reg.tif']);
    Grid_488_raw = bfGetPlane(r_488, 1);
    Grid_647_reg = bfGetPlane(r_647, 1);
    load([selpath, '/', 'Corrections', '/', 'OS_Calibration_File.mat'])
end

%Splits flat feilds and generates flat fieds correction files
if ~isfile([selpath, '/', 'Corrections', '/', 'OS_FF647.tif'])
    [avg_FFC_488, avg_FFC_647] = FFC_splitter(selpath, OS_Calibration_File);
else
    r_488 = bfGetReader([selpath, '/', 'Corrections', '/', 'OS_FF488.tif']);
    r_647 = bfGetReader([selpath, '/', 'Corrections', '/', 'OS_FF647.tif']);
    avg_FFC_488 = bfGetPlane(r_488, 1);
    avg_FFC_647 = bfGetPlane(r_647, 1);
end

%Splits live cell imaging data
if ~isfile([selpath, '/', 'Corrections', '/','Splitting_done.txt'])
    IMG_splitter(selpath, OS_Calibration_File);
end

%Performs background subtraction, FFC, and bleach correction (simple ratio)
if ~isfile([selpath, '/', 'Corrections', '/','Corrections_done.txt'])
    IMG_processing(selpath, avg_FFC_488, avg_FFC_647);
end

%Image (647 channel) registration
if ~isfile([selpath, '/', 'Corrections', '/','647_Registration_done.txt'])
    IMGREG_V2(selpath, Grid_488_raw, Grid_647_reg);
end

%Interpolation correction (488 Channel) and dz channel generation
if ~isfile([selpath, '/', 'Corrections', '/','dz_generation_done.txt'])
    dz_channel_generator(selpath, y);
end

%Organizes data to be ready for CMEanalysis
if sort4CMEanalysis == 1
    if ~isfile([selpath, '/', 'Corrections', '/','organization_done.txt'])
        data_org_for_CMEanalysis(selpath);
    end
end

delete(poolobj)
clc; fprintf('Data ready for transfer. Have a good day :)\n');

if system_shutdown == 1
    system('shutdown -s')
end
end