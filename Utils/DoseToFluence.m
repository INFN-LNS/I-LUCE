% Dose [Gy] to Fluence [cm-2] conversion for a given particle
% for a given stopping power
%

% Dose in Gy
%
Dose = input('Insert a dose value in Gy ');

% Stopping power in MeV*cm2/g
%
%StoppingPower = 95.89; %@3.9 MeV  in Water
StoppingPower = 75.59; %@5.3 MeV  in Water

% Converting the stopping power in J*cm2/Kg
%
StoppingPower = StoppingPower*1.60218e-13/1e-3;

% Fluence = Dose/SoppingPower
%
Fluence = Dose/StoppingPower;

disp('The partile fluence is ' + string(Fluence) + ' cm-2')