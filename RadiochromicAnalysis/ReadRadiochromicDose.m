% NAME: ReadRadiochromiDose
%
% Read an RCF and extract the dose in a ROI applying a given calibration
% fit already generated
%
%%
close all
clear all

%% Add folders
addpath("Functions/");

%% Main while loop
%
G = 0;

while G == 0
    disp('Have you a calibration file already prepared?')
    disp('y/n or press x to exit')

    P = input('Enter answer: ', 's');

    switch P
        case 'y'
%% Read the calibration file with Doses, net Opetical Densities and Fit parameters.
% these info comes from the procedure that read the calibration films and
% perform the polynomial fit
%
[FitParameters, Doses, NetOpticalDensities] = ReadCalibrationDataAndFitParameter();


%% Extraction of the dose value
%
% Procedure for the extraction of the dose
% 1.- Open a file with a RCF image
% 2.- Extract the average pixel values of a given ROI for the RED channel
% of the image
% 3.- Calculate the 'net Optical density' (netOD) divinding the average value in
% the ROI by the background
% 4.- Insert the netOD in the polynomial formula to extract the correct
% dose value



G = G + 1;

 case 'n'
    disp('Enter the three fit parameters in the order P1, P2 and P3')

    G = G + 1;

    case 'x'
            disp('The program has been terminated')
            G = G + 1;
        otherwise
            disp('Insert a valid value')
    end
    return
end


