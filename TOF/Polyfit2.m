%%
data = textread('/Users/Ferdi/OneDrive/Documenti/MATLAB/Prova FitPolinomiale/EneCin_EneRes_Protons_50umAl.txt','', 'delimiter', '.', ... 
                'headerlines',1,'emptyvalue', NaN);
%%
x= data(:,1);
y= data(:,2);

%xi=[0:0.1:4];
figure;
plot(x,y,'.')
hold on
%%
p = polyfit(x,y,9)
yi1=polyval(p,x)
format compact
plot(x,y,'r.',x,yi1,'b.')