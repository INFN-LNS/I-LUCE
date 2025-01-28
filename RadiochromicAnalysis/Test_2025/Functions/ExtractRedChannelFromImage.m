function [ImageRedLayer] = ExtractRedChannelFromImage(Image)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


ImageRedLayer = Image;
ImageGreenLayer = Image;
ImageBlueLayer = Image;

ImageRedLayer(:,:,2) = 0; 
ImageRedLayer(:,:,3) = 0; 
end