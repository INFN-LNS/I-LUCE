function [DoseValues, NetOpticalDensities] = OpenAndReadCalibrationCurve()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Read a file with the calibration curve already acquired
% The file to be read contains two columns.
% Column 1 is the dose
% Columns 2 is the NET Optical density already calculated
% The file must be written with this same script
%
filter = {'*.txt'};
[AcquiredCalibrationCurve, path] = uigetfile(filter,...
    'File Selection','CalibrationData.txt');

ReadCalibrationCurve = load(strcat(path,AcquiredCalibrationCurve));
DoseValues = ReadCalibrationCurve(:,1);
NetOpticalDensities = ReadCalibrationCurve(:,2);


end