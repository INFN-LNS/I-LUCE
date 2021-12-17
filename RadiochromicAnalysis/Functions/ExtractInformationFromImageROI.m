function [AveragePixelsValue,StandardDeviationPixelsValue, ...
    HorizontalDimensionInPixels, VerticalDimensionInPixels] = ExtractInformationFromImageROI(Image)

% ExtractInformationFromImageROI: Extract relevant information from an image ROI
%
% GAP Cirrone, December 2021
%
% Input: Image
% Output: AveragePixelsValue
%         StandardDeviationPixelsValue: Standard deviation of the pixels values
%         VerticalDimension: vertical dimension (in pixels) of the ROI
%         HorizontalDimension: horizontal dimension (in pixels) of the ROI

% Display the image
%
% As the image is an RGB image we have to extract only the first layer,
% corresponding to the RED
%
Image = Image(:,:,1);
Figure1 = figure('name', 'Downloaded image');
imshow(Image, 'InitialMagnification',600);

% Asks for a ROI
%
disp('Select a ROI on the picture using mouse');

ImageROI = drawrectangle('Label','ROI');

% Extracts the pixels value of the ROI

% Rounds the ROI vertexs
Vertex1 = round(ImageROI.Vertices(1,:));
Vertex2 = round(ImageROI.Vertices(2,:));
Vertex3 = round(ImageROI.Vertices(3,:));
Vertex4 = round(ImageROI.Vertices(4,:));

HorizontalDimensionInPixels = round(ImageROI.Vertices(3,1)) - round(ImageROI.Vertices(2,1));
VerticalDimensionInPixels =  round(ImageROI.Vertices(3,2)) - round(ImageROI.Vertices(1,2));



%Finally close the Figure showing the background image and ROI
%
close(Figure1)

% Calculate parameters of the extracted ROI for the Background

% Extract from the background red image the part corresponding to the ROI
%
ExtractedImage = Image(Vertex1(2):Vertex3(2),...
    Vertex1(1):Vertex3(1));

% Calculate the average of the pixels values
%
AveragePixelsValue = mean2(ExtractedImage);

% Calculate the standard deviation of the pixels values
%
StandardDeviationPixelsValue = std2(ExtractedImage);
end