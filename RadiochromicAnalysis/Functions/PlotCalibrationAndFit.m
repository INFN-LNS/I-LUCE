function PlotCalibrationAndFit(NetOpticalDensities, DosesValues, fitresult)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

xData = NetOpticalDensities;
yData = DosesValues;

% Create figure
figure1 = figure('Name','3rd polynomial fit','Color',[1 1 1], 'Position', [10 10 1100 800]);

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% Create plot with two data set: the calibration data and the fit
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
h3 = plot(xDataNew, fit_high);

% Fix only the lower limit 
%
xlim([0,inf]);
ylim([0,inf]);

% legend of the plot
%
legend1 = legend([h1 h2],{'Calibration data','3rd order polynomial fit'});
legend1.Location = 'northwest';
legend1.Interpreter = 'latex';
set(legend1, 'Box', 'on', 'Color', [0.8,0.8,0.8], 'EdgeColor', get(legend1, 'Color' )) ;
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
set(axes1,'FontSize',18, 'LineWidth', 2); 



% Add the information of the fit result inside the plot
% This shall be moved in the plot header
%
txt = strcat('The fit equation is: P3*X^3 + P2*X^2 + P1*X');
txt1 = strcat('P3 =  ', num2str(fitresult.P3));
txt2 = strcat('P2 =  ', num2str(fitresult.P2));
txt3 = strcat('P1 =  ', num2str(fitresult.P1));

text(0.05, 3, txt,'FontSize',16,'Interpreter','none')
text(0.05, 3.5, txt1,'FontSize',16,'Interpreter','latex')
text(0.05, 4.0, txt2,'FontSize',16,'Interpreter','latex')
text(0.05, 4.5, txt3,'FontSize',16, 'Interpreter','latex')

