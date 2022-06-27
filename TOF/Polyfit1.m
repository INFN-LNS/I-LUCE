%%
data = textread('/Users/Ferdi/OneDrive/Documenti/MATLAB/Prova FitPolinomiale/EneCin_Res_6.5_Al_0.5_4MeV_protons.txt','', 'delimiter', '.', ... 
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