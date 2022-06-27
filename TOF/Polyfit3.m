%%
data = textread('/Users/Ferdi/OneDrive/Documenti/MATLAB/Prova FitPolinomiale/EneRes_EneLoss_100umDD.txt','', 'delimiter', '.', ... 
                'headerlines',1,'emptyvalue', NaN);
%%
x= data(:,1);
y= data(:,2);

%xi=[0:0.1:4];
figure;
plot(x,y,'o')
xlim([4 22])
hold on
%%
p = polyfit(x,y,9)
yi1=polyval(p,x)
format compact
plot(x,y,'ro',x,yi1,'b.')