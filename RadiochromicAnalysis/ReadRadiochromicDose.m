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

            %% Read the calibration file with Doses, net Optical Densities and Fit parameters.
            % these info comes from the procedure that read the calibration films and
            % perform the polynomial fit
            %y
            filter = {'*.txt'};
            [FileNameCalibrationDataAndFit, path] = uigetfile(filter,...
                    'File Selection','CalibrationDataAndFit.txt');

                % Open the file just called
                %
                pathForCalibrationDataAndFit = strcat(path, FileNameCalibrationDataAndFit);

                [P1, P2, P3, BackgroundValue,...
                    Doses,...
                    NetOpticalDensities] = ReadCalibrationDataAndFit(pathForCalibrationDataAndFit);
            %% Measure of the background
            %
            disp('Open one RCF image to measure the background.')
            [filename,user_canceled] = imgetfile;

            % Read the image
            %
            Image = imread(filename);
            ImageRed = ExtractRedChannelFromImage(Image);

            [MeasuredAverageBackground, StandardDeviationPixelsValueBackground] = ...
                ExtractInformationFromImageROI(ImageRed);

            disp('The measured average background is: ' + string(MeasuredAverageBackground));


            %% Open the image to be analised
            %

            disp('Open one RCF image.')
            [filename,user_canceled] = imgetfile;

            % Read the image
            %
            Image = imread(filename);
            ImageRed = ExtractRedChannelFromImage(Image);

            [AveragePixelsValue, ...
                StandardDeviationPixelsValue,...
                HorizontalDimension,...
                VerticalDimension] = ...
                ExtractInformationFromImageROI(ImageRed);

            %% Calculation of the Net Optical Density 
            %
            NetOpticalDensity = -log10(AveragePixelsValue/MeasuredAverageBackground);

            %% Extraction of the dose value using the fit parameters
            %
            Dose = P3*NetOpticalDensity.^3 + ...
                P2*NetOpticalDensity.^2 + ...
                P1*NetOpticalDensity;

             disp('The stimated dose is ' + string(Dose) + ' cGy')

     
            G = G + 1;

        case 'n'
            disp('Enter the three fit parameters in the order P1, P2 and P3')
            FitParameterAsInput = [];
            for i = 1:1:3
                FitParameter = input('Insert the fit parameter P' + string(i) + ':');
                FitParametersAsInput(i) = FitParameter;
            end

            %% Measure of the background
            %
            disp('Open one RCF image to measure the background.')
            [filename,user_canceled] = imgetfile;

            % Read the image
            %
            Image = imread(filename);
            ImageRed = ExtractRedChannelFromImage(Image);

            [MeasuredAverageBackground, StandardDeviationPixelsValueBackground] = ...
                ExtractInformationFromImageROI(ImageRed);


 %           BackgroundValueAsInput = input('Now eneter the measured value of the background: ');


             %% Open the image to be analised
             %
            disp('Open one RCF image.')
            [filename,user_canceled] = imgetfile;

            % Read the image
            %
            Image = imread(filename);
            ImageRed = ExtractRedChannelFromImage(Image);

            [AveragePixelsValue, ...
                StandardDeviationPixelsValue,...
                HorizontalDimensionInPixels,...
                VerticalDimensionInPixels] = ...
                ExtractInformationFromImageROI(ImageRed);

            %% Calculation of the Net Optical Density 
            %
            NetOpticalDensity = -log10(AveragePixelsValue/MeasuredAverageBackground);

            %% Extraction of the dose (in cGy) value using the fit parameters
            %
            Dose = FitParametersAsInput(3)*NetOpticalDensity.^3 + ...
                FitParametersAsInput(2)*NetOpticalDensity.^2 + ...
                FitParametersAsInput(1)*NetOpticalDensity;

            disp('The stimated dose is ' + string(Dose) + ' cGy')

            % DPI (dots per inch) of the scanned image
            %
            DPI = 600;
            % Pixel dimension in cm
            %
            pixelDimension = 2.54/600;
            
            % Surface of the ROI over which the dose was calculated in cm2
            %
            ROISurface = HorizontalDimensionInPixels*pixelDimension*VerticalDimensionInPixels*pixelDimension;

            disp('ROI surface over which the dose is calculated: ' + string(ROISurface) + ' cm2')

            
            G = G + 1;

        case 'x'
            disp('The program has been terminated')
            G = G + 1;
        otherwise
            disp('Insert a valid value')
    end
    return
end


