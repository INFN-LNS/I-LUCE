function [filenameCalibrationCurveAndFit] = WriteDataAndFitParameters(DosesValues, NetOpticalDensities, fitresult)
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

% Write the file with the calibration points and FIT results
%
fprintf(fid,'%6s %12s\n','Dose [ Gy ]','Net OD');
fprintf(fid, '%6.2f %12.8f\n',[DosesValues; NetOpticalDensities]);
fprintf(fid, '%12s\n', '------------------------ ');
fprintf(fid, '%12s\n', 'POLYNOMIAL FIT PARAMETERS ');
fprintf(fid, '%12s\n', '                   ');
fprintf(fid, '%12s\n', 'The FIT formula is: P3*X^3 + P2*X^2 + P1*X');
fprintf(fid, '%12s\n', 'Values of the fit parameters:');
fprintf(fid, '%5s %5s\n', 'P3= ', num2str(fitresult.P3));
fprintf(fid, '%5s %5s\n', 'P2= ', num2str(fitresult.P2));
fprintf(fid, '%5s %5s\n', 'P1= ', num2str(fitresult.P1));
fclose(fid);
end

