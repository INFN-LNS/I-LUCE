function [h] = PlotCalibrationCurve(NetOpticalDensities, DosesValues)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

FigureDoseResponseCurve = figure('Color',[1 1 1]);

% Create axes
axes1 = axes('Parent',FigureDoseResponseCurve);
hold(axes1,'on');

h = plot(NetOpticalDensities, DosesValues,'MarkerSize',16,...
    'Marker','o','LineStyle','none','LineWidth', 4);

xlabel('Net Optical density')
ylabel('Dose [Gy]')

box(axes1,'on');
hold(axes1,'off');

% Set the remaining axes properties
set(axes1,'FontSize',18,'GridColor',[0.8 0.8 0.8],'LineWidth',2,'XGrid',...
    'on','YGrid','on');
end