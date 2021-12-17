function [FitParameters, Doses, NetOpticalDensities, BackgroundValue] = ReadCalibrationDataAndFitParameter()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


filter = {'*.txt'};
[FileNameCalibrationAndFitParameters, path] = uigetfile(filter,...
    'File Selection','CalibrationDataAndFit.txt');

% Open the file just called
%
pathForCalibrationAndFitParameters = strcat(path, FileNameCalibrationAndFitParameters);

% open the file
%
fid=fopen(pathForCalibrationAndFitParameters); 

FitParameters = textscan(fid,'%f%f%f', 'Delimiter','','headerlines',4);  

BackgroundValue = textscan(fid,'%12.f', 'headerlines',1);
DosesAndONetpticalDensities = textscan(fid,'%12.f%12.f', 'Delimiter','', 'headerlines',1);


% Extraction from the cells the vector with the elevant data
%
FitParameters = [FitParameters{1} FitParameters{3} FitParameters{3}]';
Doses = DosesAndONetpticalDensities{1}(:);
NetOpticalDensities = DosesAndONetpticalDensities{2}(:);
BackgroundValue = BackgroundValue{1};

fclose(fid);

end