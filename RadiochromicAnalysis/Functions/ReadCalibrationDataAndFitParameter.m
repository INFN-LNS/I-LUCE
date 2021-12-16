function [FitParameters, Doses, NetOpticalDensities] = ReadCalibrationDataAndFitParameter()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


filter = {'*.txt'};
[FileNameCalibrationAndFitParameters, path] = uigetfile(filter,...
    'File Selection','CalibrationDataAndFit.txt');

% Open the file just called
%
pathForCalibrationAndFitParameters = strcat(path, FileNameCalibrationAndFitParameters);
% open the file
fid=fopen(pathForCalibrationAndFitParameters); 

line = 5;
FitParameters = textscan(fid,'%f %f %f', 1, 'Delimiter', '', 'WhiteSpace', '', 'headerlines',line-1);  % Read 2 header lines
DosesAndONetpticalDensities = textscan(fid,'%f%f', 'Delimiter',{'','\n'}, 'headerlines',3);

% Extraction from the cells the vector with the elevant data
%
FitParameters = [FitParameters{1} FitParameters{3} FitParameters{3}]';
Doses = DosesAndONetpticalDensities{1}(:);
NetOpticalDensities = DosesAndONetpticalDensities{2}(:);

fclose(fid);

end