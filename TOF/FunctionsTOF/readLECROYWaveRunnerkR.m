function [tabularData] = readLECROYWaveRunnerkR(filename)
%% Import data from text file
% Script for importing data from the following text file:
%
%    filename: /home/pablo.cirrone/WaterTarget/2022-05-09/Manual/C1--XX--00000.txt
%
% Auto-generated by MATLAB on 10-May-2022 16:08:26

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 2);

% Specify range and delimiter
opts.DataLines = [6, Inf];
opts.Delimiter = " ";

% Specify column names and types
opts.VariableNames = ["Time", "Ampl"];
opts.VariableTypes = ["double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts.ConsecutiveDelimitersRule = "join";
opts.LeadingDelimitersRule = "ignore";

% Import the data

tabularData = readtable("filename", opts);
%tabularData = readtable("/home/pablo.cirrone/WaterTarget/2022-05-09/Manual/C1--XX--00000.txt", opts);


%% Clear temporary variables
clear opts