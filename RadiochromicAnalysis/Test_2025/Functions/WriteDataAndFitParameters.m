function [filenameCalibrationCurveAndFit] = ...
    WriteDataAndFitParameters(DoseValues,...
    NetOpticalDensities,...
    fitresult,...
    AveragePixelsValueBackground)


%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

% Write calibration points and fit parameters in the same file
%
disp(['Choose where the calibration data (Dose and NetOD)...' ...
    ' and Fit results will be saved'])

% Open a dialog box specifying the file name where to save the calibration
% data
% A default file name is given
%
filter = {'*.txt'};
[filenameCalibrationCurveAndFit,path] = uiputfile(filter,...
    'File Selection','CalibrationDataAndFit.txt');

% Open the file just defined
%
pathForCalibrationCurveAndFitParameters= strcat(path, filenameCalibrationCurveAndFit);

fid = fopen(pathForCalibrationCurveAndFitParameters, 'wt');

fprintf(fid, '%12s\n', 'POLYNOMIAL FIT PARAMETERS ');
fprintf(fid, '%12s\n', 'The FIT formula is: P3*X^3 + P2*X^2 + P1*X');
fprintf(fid, '\n');
fprintf(fid,'%6s %6s %6s\n','P1','P2','P3');
fprintf(fid, '%f %f %f\n',[fitresult.P1, fitresult.P2, fitresult.P3]);
fprintf(fid, '\n');
fprintf(fid, '%16s\n', 'Background value');
fprintf(fid, '%6f\n', AveragePixelsValueBackground);
fprintf(fid, '\n');
fprintf(fid,'%6s %6s\n','Dose [ Gy ]','Net OD');
fprintf(fid, '%f %f\n',[DoseValues, NetOpticalDensities]'); 

fclose(fid);
end

