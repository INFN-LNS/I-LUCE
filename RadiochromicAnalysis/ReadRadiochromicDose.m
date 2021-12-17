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
            disp('Open the file with the calibration data and fit parameters.')

            %% Read the calibration file with Doses, net Opetical Densities and Fit parameters.
            % these info comes from the procedure that read the calibration films and
            % perform the polynomial fit
            %
            [FitParameters, Doses, NetOpticalDensities, BackgroundValue] = ReadCalibrationDataAndFitParameter();

            %% Open the image to be analised
            %

            disp('Open one RCF image.')
            [filename,user_canceled] = imgetfile;

            % Read the image
            %
            Image = imread(filename);
            ImageRed = ExtractRedChannelFromImage(Image);

            [AveragePixelsValue, StandardDeviationPixelsValue] = ...
                ExtractInformationFromImageROI(ImageRed);

            %% Calculation of the Net Optical Density 
            %
            NetOpticalDensity = -log10(AveragePixelsValue/BackgroundValue);

            %% Extraction of the dose value using the fit parameters
            %
            Dose = FitParameters(3)*NetOpticalDensity.^3 + ...
                FitParameters(2)*NetOpticalDensities.^2 + ...
                FitParameters(1)*NetOpticalDensities;
             
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
            FitParameterAsInput = [];
            for i = 1:1:3
                FitParameter = input('Insert the fit parameter P' + string(i) + ':');
                FitParametersAsInput(i) = FitParameter;
            end

            BackgroundValueAsInput = input('Now eneter the measured value of the background: ');


             %% Open the image to be analised
             %

            disp('Open one RCF image.')
            [filename,user_canceled] = imgetfile;

            % Read the image
            %
            Image = imread(filename);
            ImageRed = ExtractRedChannelFromImage(Image);

            [AveragePixelsValue, StandardDeviationPixelsValue] = ...
                ExtractInformationFromImageROI(ImageRed);

            %% Calculation of the Net Optical Density 
            %
            NetOpticalDensity = -log10(AveragePixelsValue/BackgroundValueAsInput);

            %% Extraction of the dose value using the fit parameters
            %
            Dose = FitParametersAsInput(3)*NetOpticalDensity.^3 + ...
                FitParametersAsInput(2)*NetOpticalDensity.^2 + ...
                FitParametersAsInput(1)*NetOpticalDensity;

            disp('The stimated dose is ' + string(Dose) + ' CGy')


            G = G + 1;

        case 'x'
            disp('The program has been terminated')
            G = G + 1;
        otherwise
            disp('Insert a valid value')
    end
    return
end


