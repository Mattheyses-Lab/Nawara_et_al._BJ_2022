%%% Function created by Nawara T. (Mattheyses lab - 10/06/2022) compatible
%%% with MATLAB R2020b
%%% Dream STAR
%%% Data spliting, corection, dz_channel generator, and data export
%%% package
%%% Suports data aquired with optosplit or dual-camera systems

%This work is licensed under the Creative Commons Attribution 4.0
%International License. To view a copy of this license, visit
%http://creativecommons.org/licenses/by/4.0/ or send a letter to Creative
%Commons, PO Box 1866, Mountain View, CA 94042, USA.

function DrSTAR

close all

path2code = pwd;

if ismac
    addpath(genpath([path2code, '/bfmatlab']));
    addpath(genpath([path2code, '/image registration']));
elseif ispc
    addpath(genpath([path2code, '\bfmatlab']));
    addpath(genpath([path2code, '\image registration']));
else
    disp('Platform not supported')
end
load('Colormaps.mat')
assignin('base','Colormaps', Colormaps)

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

% Specify whetehr data is aquired using Optosplit (OS) or dual-camera system (DC)
if ~isfile([selpath, '/', 'Corrections', '/', 'data_type.mat'])
    data_type = input('Data aquierd using Optosplit or Dual-camera system (OS/DC)?   ','s');
    data_type_path = [selpath, '/', 'Corrections', '/', 'data_type.mat'];
    save(data_type_path,'data_type');
else
    load([selpath, '/', 'Corrections', '/', 'data_type.mat'])
end

switch data_type
    case 'OS'
        %Grid splitting and Optosplit calibration
        if ~isfile([selpath, '/', 'Corrections', '/', 'Calibration_File.mat'])
            [Calibration_File, Grid_488_raw, Grid_647_reg] = OS_Cali(selpath);
        else
            r_488 = bfGetReader([selpath, '/', 'Corrections', '/', 'Grid_488_raw.tif']);
            r_647 = bfGetReader([selpath, '/', 'Corrections', '/', 'Grid_647_reg.tif']);
            Grid_488_raw = bfGetPlane(r_488, 1);
            Grid_647_reg = bfGetPlane(r_647, 1);
            load([selpath, '/', 'Corrections', '/', 'Calibration_File.mat'])
        end
        
    case 'DC'
        %Grid splitting and Optosplit calibration
        if ~isfile([selpath, '/', 'Corrections', '/', 'Calibration_File.mat'])
            [Calibration_File, Grid_488_raw, Grid_647_reg] = DC_Cali(selpath);
        else
            r_488 = bfGetReader([selpath, '/', 'Corrections', '/', 'Grid_488_raw.tif']);
            r_647 = bfGetReader([selpath, '/', 'Corrections', '/', 'Grid_647_reg.tif']);
            Grid_488_raw = bfGetPlane(r_488, 1);
            Grid_647_reg = bfGetPlane(r_647, 1);
            load([selpath, '/', 'Corrections', '/', 'Calibration_File.mat'])
        end
end
            

        
%Splits flat feilds and generates flat fieds correction files
if ~isfile([selpath, '/', 'Corrections', '/', 'FF647.tif'])
    [avg_FFC_488, avg_FFC_647] = FFC_splitter(selpath, Calibration_File);
else
    r_488 = bfGetReader([selpath, '/', 'Corrections', '/', 'FF488.tif']);
    r_647 = bfGetReader([selpath, '/', 'Corrections', '/', 'FF647.tif']);
    avg_FFC_488 = bfGetPlane(r_488, 1);
    avg_FFC_647 = bfGetPlane(r_647, 1);
end
        
switch data_type
    case 'OS'
        %Splits live cell imaging data
        if ~isfile([selpath, '/', 'Corrections', '/','Splitting_done.txt'])
            OS_IMG_splitter(selpath, Calibration_File);
        end
        
    case 'DC'
        %Splits live cell imaging data
        if ~isfile([selpath, '/', 'Corrections', '/','Splitting_done.txt'])
            DC_IMG_splitter(selpath, Calibration_File);
        end
end

%Performs background subtraction, FFC, and bleach correction (simple ratio)
if ~isfile([selpath, '/', 'Corrections', '/','Corrections_done.txt'])
    IMG_processing(selpath, avg_FFC_488, avg_FFC_647);
end

%Image (647 channel) registration
if ~isfile([selpath, '/', 'Corrections', '/','647_Registration_done.txt'])
    IMGREG_V2(selpath, Grid_488_raw, Grid_647_reg, Colormaps);
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
