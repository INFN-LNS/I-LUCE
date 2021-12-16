% GAP Cirrone, November 2021
% Last revision: December 2021
%
% Calibration in dose of radiochromic films
%

%% Close and clean all variables and pictures.
%
clear all
close all

%% Add folders
addpath("Functions/");

%%  main while loop
%
G = 0;

while G == 0
    disp('Have you a calibration already acquired?')
    disp('y/n or press x to exit')

    P = input('Enter answer: ', 's');

    switch P
        case 'y'
            % Open and read an old calibration file
            %
            [DoseValues,NetOpticalDensities] = OpenAndReadCalibrationCurve();

            % PlotCalibrationCurve(NetOpticalDensities, DosesValues) and
            % makes a 3rd order polynomial fit
            %
            [fitresult, gof] = RCFPolynomialFit(NetOpticalDensities, DoseValues);

            % Plot the calibration data toegheter with the polynomial fit
            % curve and make a separate plot with the residuals
            %
            PlotCalibrationAndFit(NetOpticalDensities, DoseValues, fitresult)

            % Write the calibration data and the results from the
            % polynomial fit in a dingle text file
            %
            WriteDataAndFitParameters(DoseValues, NetOpticalDensities', fitresult)

            G = G + 1;

        case 'n'
            % Enter in the case 'n' a previous calibration 
            % is not present and data must be acquired 
            % reading a set of irradiated RCF
            %
            disp('Have you a file containing a list of calibration doses?')
            disp('y/n or press x to exit')

            Q = input('Enter answer: ', 's');
            if Q == 'y'
                % Here a file with the calibration dose values can be
                % uploaded
                %
                disp('Please, load a file with the dose values');
                filter = {'*.txt'};
                [filenameDoseValuesForCalibration, path] = uigetfile(filter,...
                    'File Selection','DoseValuesForCalibration.txt');

                % Open the file just called
                %
                pathForDoseValuesForCalibration = strcat(path, filenameDoseValuesForCalibration);

                % Read the files with the doses
                %
                DoseValues = load(pathForDoseValuesForCalibration);
                DosesPoints = length(DoseValues);
                
                %Load the background image and extract the RED channel
                %
                % Asks for opening an image
                %
                disp('Open the Background image')
                [filename,user_canceled] = imgetfile;

                % Read the image
                %
                BackgroundImage = imread(filename);

                [BackgroundImageRedLayer] = ExtractRedChannelFromImage(BackgroundImage);

                % Asks for a folder where save the red image of the background
                %
                disp('Choose where the red image has to be saved')

                filter = {'*.TIFF'};
                [filenameImage_red,path] = uiputfile(filter,...
                    'File Selection','Fondo.TIFF');

                % Save the Background image in the red channel
                %
                pathToSave = convertCharsToStrings(convertCharsToStrings(path) + ...
                    filenameImage_red);
                imwrite(BackgroundImageRedLayer, pathToSave);
                disp('The image has been saved in ' + convertCharsToStrings(filenameImage_red));

                %% Open the red image, ROI it and extract information
                %
                disp('Open the Background red image')
                [filename,user_canceled] = imgetfile;

                % Read the background image
                BackgroundImage = imread(filename);

                % The following funtion extracts the pixel values
                % and the correposnding standard deviation from a ROI of a given image
                %
                [AveragePixelsValueBackground, StandardDeviationPixelsValueBackground] = ...
                    ExtractInformationFromImageROI(BackgroundImage);

                %% Asks to open the calibration images and stores the ROI information
                %
                AveragePixelsValueCalibration = [];
                StandardDeviationPixelsValueCalibration = [];

                for i = 1 : 1 : DosesPoints
                    disp('Open the RCF corresponding to the Dose ' +  string(i) + ': ')
                    [filename,user_canceled] = imgetfile;

                    % Read the image
                    %
                    Image = imread(filename);
                    ImageBackground = ExtractRedChannelFromImage(Image);

                    [AveragePixelsValue, StandardDeviationPixelsValue] = ...
                        ExtractInformationFromImageROI(ImageBackground);

                    AveragePixelsValueCalibration(i) = AveragePixelsValue;
                    StandardDeviationPixelsValueCalibration(i) = StandardDeviationPixelsValue;
                end

                %% Calculate the Net Optical Density
                %
                NetOpticalDensities = -log10(AveragePixelsValueCalibration/AveragePixelsValueBackground);
                %% Write in a file the Calibration points
                % (Doses and the NetOpticalDensities)
                %
                disp('Choose where the calibration data (Dose and NetOD) will be saved')

                % Open a dialog box specifying the file name where to save the calibration
                % data
                % A default file name is given
                %

                % Asks to open a file where to write the calibration points
                %
                filter = {'*.txt'};
                [filenameCalibrationCurve, path] = uiputfile(filter,...
                    'File Selection','CalibrationData.txt');

                pathForCalibrationFile = strcat(path, filenameCalibrationCurve);

                % Open the file just defined
                %
                fid = fopen(pathForCalibrationFile, 'wt');

                % Write the file with the calibration points
                %
                fprintf(fid, '%f %f\n', [DoseValues'; NetOpticalDensities]);
                fclose(fid);

                %% Plot the calibration curve
                %
                %PlotCalibrationCurve(NetOpticalDensities, DosesValues)

                % PlotCalibrationCurve(NetOpticalDensities, DosesValues) and
                % makes a 3rd order polynomial fit
                %
                [fitresult, gof] = RCFPolynomialFit(NetOpticalDensities, DoseValues);

                % Plot the calibration data toegheter with the polynomial fit
                % curve and make a separate plot with the residuals
                %
                PlotCalibrationAndFit(NetOpticalDensities, DoseValues', fitresult)

                % Write the calibration data and the results from the
                % polynomial fit in a dingle text file
                %
                WriteDataAndFitParameters(DoseValues, NetOpticalDensities', fitresult, AveragePixelsValueBackground)
            end
                
               
            if Q == 'n'
                % In this case, the file with the used doses is not
                % available and the User have to explicitly write down the
                % values of the doses used for the calibration
                %
                %% Asks for the number of dose points
                %
                prompt = 'Enter the number of dose points:  ';
                DosesPoints = input(prompt);

                %% Enter the dose values for each point
                disp('Now, enter the doses values in Gy (press return after each value)')

                DosesValues = [];
                for i =  1 : 1 : DosesPoints
                    prompt = 'Enter the dose value number  ' +  string(i) + ': ';
                    DosesValues(i) = input(prompt);
                end

                %% Load the background image and extract the RED channel
                %
                % Asks for opening an image
                %
                disp('Open the Background image')
                [filename,user_canceled] = imgetfile;

                % Read the image
                %
                BackgroundImage = imread(filename);

                [BackgroundImageRedLayer] = ExtractRedChannelFromImage(BackgroundImage);

                % Asks for a folder where save the red image of the background
                %
                disp('Choose where the red image has to be saved')

                filter = {'*.TIFF'};
                [filenameImage_red,path] = uiputfile(filter,...
                    'File Selection','Fondo.TIFF');

                % Save the Background image in the red channel
                %
                pathToSave = convertCharsToStrings(convertCharsToStrings(path) + ...
                    filenameImage_red);
                imwrite(BackgroundImageRedLayer, pathToSave);
                disp('The image has been saved in ' + convertCharsToStrings(filenameImage_red));

                %% Open the red image, ROI it and extract information
                %
                disp('Open the Background red image')
                [filename,user_canceled] = imgetfile;

                % Read the background image
                BackgroundImage = imread(filename);

                % The following funtion extracts the pixel values
                % and the correposnding standard deviation from a ROI of a given image
                %
                [AveragePixelsValueBackground, StandardDeviationPixelsValueBackground] = ...
                    ExtractInformationFromImageROI(BackgroundImage);

                %% Asks to open the calibration images and stores the ROI information
                %
                AveragePixelsValueCalibration = [];
                StandardDeviationPixelsValueCalibration = [];

                for i = 1 : 1 : DosesPoints
                    disp('Open the RCF corresponding to the Dose ' +  string(i) + ': ')
                    [filename,user_canceled] = imgetfile;

                    % Read the image
                    %
                    Image = imread(filename);
                    ImageBackground = ExtractRedChannelFromImage(Image);

                    [AveragePixelsValue, StandardDeviationPixelsValue] = ...
                        ExtractInformationFromImageROI(ImageBackground);

                    AveragePixelsValueCalibration(i) = AveragePixelsValue;
                    StandardDeviationPixelsValueCalibration(i) = StandardDeviationPixelsValue;
                end

                %% Calculate the Net Optical Density
                %
                NetOpticalDensities = -log10(AveragePixelsValueCalibration/AveragePixelsValueBackground);
                %% Write in a file the Calibration points
                % (Doses and the NetOpticalDensities)
                %
                disp('Choose where the calibration data (Dose and NetOD) will be saved')

                % Open a dialog box specifying the file name where to save the calibration
                % data
                % A default file name is given
                %

                % Asks to open a file where to write the calibration points
                %
                filter = {'*.txt'};
                [filenameCalibrationCurve, path] = uiputfile(filter,...
                    'File Selection','CalibrationData.txt');

                pathForCalibrationFile = strcat(path, filenameCalibrationCurve);

                % Open the file just defined
                %
                fid = fopen(pathForCalibrationFile, 'wt');

                % Write the file with the calibration points
                %
                fprintf(fid, '%f %f\n', [DosesValues; NetOpticalDensities]);
                fclose(fid);

                %% Plot the calibration curve
                %
                %PlotCalibrationCurve(NetOpticalDensities, DosesValues)

                % PlotCalibrationCurve(NetOpticalDensities, DosesValues) and
                % makes a 3rd order polynomial fit
                %
                [fitresult, gof] = RCFPolynomialFit(NetOpticalDensities, DosesValues);

                % Plot the calibration data toegheter with the polynomial fit
                % curve and make a separate plot with the residuals
                %
                PlotCalibrationAndFit(NetOpticalDensities, DosesValues, fitresult)

                % Write the calibration data and the results from the
                % polynomial fit in a dingle text file
                %
                WriteDataAndFitParameters(DosesValues, NetOpticalDensities, fitresult)
            end
        case 'x'
            disp('The program has been terminated')
            G = G + 1;
        otherwise
            disp('Insert a valid value')
    end
    return
end




