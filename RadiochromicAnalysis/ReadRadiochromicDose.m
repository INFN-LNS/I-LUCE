% NAME: ReadRadiochromiDose
%
% Read an RCF and extract the dose in a ROI applying a given calibration
% fit already generated
%
%% Read the calibration file with fitting 
%
% Function: ReadCalibrationAndFit
%

DosePoints = 6;

filter = {'*.txt'};
[FileNameCalibrationAndFitParameters, path] = uigetfile(filter,...
    'File Selection','CalibrationDataAndFit.txt');

% Open the file just called
%

pathForCalibrationAndFitParameters = strcat(path, FileNameCalibrationAndFitParameters);

opts = detectImportOptions(pathForCalibrationAndFitParameters);
opts.DataLines = [2 7];
opts.VariableNames = {'Dose','Net OD'};
opts.Delimiter(' ')
T_first = readtable(pathForCalibrationAndFitParameters,opts)

