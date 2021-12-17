% Fluence [cm-2] from dose [Gy]
%

% Dose in Gy
%
Dose = input('Insert a dose value in Gy ');

% Stopping power in MeV*cm2/g
%
StoppingPower = 146.5;

% Converting the stopping power in J*cm2/Kg
%
StoppingPower = StoppingPower*1.60218e-13/1e-3;

% Fluence = Dose/SoppingPower
%
Fluence = Dose/StoppingPower;

disp('The partile fluence is ' + string(Fluence) + ' cm-2')