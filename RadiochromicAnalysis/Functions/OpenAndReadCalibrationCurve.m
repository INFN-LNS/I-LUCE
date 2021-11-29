function [DosesValues,NetOpticalDensities] = OpenAndReadCalibrationCurve()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Read a file with the calibration curve already acquired
% The file to be read contains two columns.
% Column 1 is the dose
% Columns 2 is the NET Optical density already calculated
% The file must be written with this same script
filter = {'*.txt'};
[AcquiredCalibrationCurve, user_canceled] = uigetfile(filter,...
    'File Selection','CalibrationData.txt');
ReadCalibrationCurve = load(convertCharsToStrings(AcquiredCalibrationCurve));
DosesValues = ReadCalibrationCurve(:,1);
NetOpticalDensities = ReadCalibrationCurve(:,2);


end