function PlotCalibrationAndFit(NetOpticalDensities, DosesValues, fitresult)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

xData = NetOpticalDensities;
yData = DosesValues;

% Create figure
figure1 = figure('Name','3rd order polynomial fit','Color',[1 1 1], 'Position', [10 10 1100 800]);

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% Create plot with two data set: the calibration data and the fit
%
% Plot the experimental data
%
h1 = plot(xData,yData);
set(h1, 'LineWidth',5,...
     'MarkerSize',16,...
     'Marker','o',...
     'LineStyle','none',...
     'Color',[0 0 1]);

% Plot the fit curve
%
h2 = plot(fitresult);
set(h2, 'LineWidth',2)

% Plot the confidence bounds 
%
% Extract the limits of the calculated coefficients
%
FitCoefficientBounds = confint(fitresult);
P1_low = FitCoefficientBounds(1);
P1_high = FitCoefficientBounds(2);
P2_low = FitCoefficientBounds(3);
P2_high = FitCoefficientBounds(4);
P3_low = FitCoefficientBounds(5);
P3_high = FitCoefficientBounds(6);

% Create a X Data series
% 
xDataNew = linspace(min(xData), max(xData));

fit_high = P3_high*xDataNew.^3 + P2_high*xDataNew.^2 + P1_high*xDataNew;
fit_low = P3_low*xDataNew.^3 + P2_low*xDataNew.^2 + P1_low*xDataNew;

h3 = plot(xDataNew, fit_high, 'b');
h4 = plot(xDataNew, fit_low, 'b');
set(h3, 'LineWidth',2);
set(h4, 'LineWidth',2);

x2 = [xDataNew, fliplr(xDataNew)];
inBetween = [fit_low, fliplr(fit_high)];

fill(x2, inBetween, 'black','FaceAlpha', '0.1');


% legend of the plot
%
legend1 = legend([h1 h2 h3 h4],{'Calibration data',...
    '3rd order polynomial fit',...
    '95$\%$ confidence bound', '95$\%$ confidence bound'});
legend1.Location = 'northwest';
legend1.Interpreter = 'latex';
set(legend1, 'Box', 'on', 'Color', [0.8,0.8,0.8]);
set(legend1, 'EdgeColor', get(legend1, 'Color' ));

%legend boxoff

box(axes1,'on');
grid(axes1,'on');
hold(axes1,'off');

% Create ylabel
ylabel('Dose [ Gy ]', 'Interpreter','latex');

% Create xlabel
xlabel('Net Optical Density [ $-\log \frac{I}{I_0}$ ]','Interpreter','latex');

set(gca,'TickLabelInterpreter','latex')

box(axes1,'on');
grid(axes1,'on');
hold(axes1,'off');

% Set the remaining axes properties
%
set(axes1,'FontSize',12, 'LineWidth', 2); 

% Add the information of the fit result in the plot title
% This shall be moved in the plot header
%
txt = strcat('Fit equation: P3*X^3 + P2*X^2 + P1*X;');
txt1 = strcat('P3 =  ', num2str(fitresult.P3));
txt2 = strcat('P2 =  ', num2str(fitresult.P2));
txt3 = strcat('P1 =  ', num2str(fitresult.P1));

% Add the plot title with information on Fit
%
title(strcat(txt, txt1,'; ',txt2,'; ',txt3));

