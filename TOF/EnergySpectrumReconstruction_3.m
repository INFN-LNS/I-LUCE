%% Experimental data 
% Non c'è un filtro e i protoni dopo il filtro si fermano dentro al rivelatore

% Rivelatore: SiC (densità=3.21 g/cm3)
% Spessore Rivelatore:10 um
% Diametro Rivelatore= 1 mm
% Distanza target-rivelatore= 20 cm

data = textread('/Users/Ferdi/OneDrive/Documenti/MATLAB/DatiProva Spettro/SenzaFiltro.txt', '', 'delimiter', '.', ... 
                'headerlines',5,'emptyvalue', NaN);

% Commenta SiC o Diamante e controlla le aree S_d 
           S='SiC';
           %S='Diamond'; %%-> display tipo di rivelatore
           dB=1; %%-> Attenuation 20 DB->10; 10 DB->3; 6 DB->2; 13 DB-> 4.45
            
  
if  strcmp(S,'SiC')
   
figure;
plot(data(:,1)*10^9,-smooth(data(:,4)*dB,10));
xlim([-100 300])
Eg=7.78*10^-6 %% [MeV] electron-hole pair energy creation 13 eV for Diamond, 7.78 for SiC 
S_d=(3.14*0.5*0.5)*10^-2; %% [cm^2]-> Display area rivelatore

else if strcmp(S,'Diamond')     
                    
figure;
plot(data(:,1)*10^9,smooth(data(:,2)*dB,10));
Eg=13*10^-6 %% [MeV] electron-hole pair energy creation 13 eV for Diamond, 7.78 for SiC 
S_d=2.5*10^-2;
end
end

%% Insert the time shift and the other parameters
dt=7.3; % [ns]

%%FLIGHT PATH%%
r=0.2; % [m]

c2=9*10^(16); % [m2/s2] c^2
c=3*10^8; % [m/s]
k=8.6*10^-5; % [eV] Cost Boltzmann
%%Type of particle
mp= 1*((938.27)/(9*10^16)); % [MeV/c2] col 4 sono alfa

e=1.6*10^-19; 

%%Flight path in cm
r_cm=r*100;

%detector solid angle
Ang=S_d/(r_cm^2); 

Dr=0.01; %% uncertainty on the distance in cm
Dt=10^-(10); % uncertainty on time in s
dAng=2*S_d*(Dr)/(r_cm)^3

tl=r/c;% time of the light in sec


 D_1=data(:,1)-(dt*10^-9)+tl; %%total time shift [sec]
 D_2=(data(:,4)*dB)
 
% Select the "time" range of interest from TOF signal using the brush tool and call the variable "proton"
 figure;
 plot(D_1*10^9,-smooth(D_2,10));
 xlim([-100 300])


%% For particle which (after the filter) stop inside the detector active layer
%%select in the plot and call the variable proton
time_min_stop=proton(1,1); % ns
Ene_time_min_stop=(0.5*mp*r^2)/((time_min_stop*10^-9)^2)

time_max_stop=102;%ns
Ene_time_max_stop=(0.5*mp*r^2)/((time_min_stop*10^-9)^2)  

D0_2(:,1)=proton(:,1);
D0_2(:,2)=proton(:,2);

t_stop=D0_2(D0_2(:,1)>time_min_stop&D0_2(:,1)<time_max_stop,1)*10^-9;
V_stop=D0_2(D0_2(:,1)>time_min_stop&D0_2(:,1)<time_max_stop,2);
figure;
plot(t_stop,V_stop);
R=50;

Dr=1;
for(i=1:length(t_stop))
 beta_2(i)=(r/(t_stop(i)))/c;
 gamma_2(i)=1/(sqrt(1-beta_2(i).^2));
 E_inc_stop(i)=(gamma_2(i)-1)*mp*c^2;
%E_stop(i)=((mp*r*r)/(2*t_stop_2(i)*t_stop_2(i)));
sigma_inc2(i)=sqrt(((mp*r/t_stop(i))^2*(Dr^2))+((((mp*r^2)/(t_stop(i).^3)).^2)*(Dt)^2));
end;
figure
plot(E_inc_stop,V_stop)

%%

for(i=1:length(t_stop))
dN_dE_stop(i)=((Eg*V_stop(i))/(e*R*E_inc_stop(i)*E_inc_stop(i)))*(0.5*t_stop(i)-Dt)
sigma_dN_stop(i)=(((Eg*V_stop(i))/(e*R*Ang))*sqrt((((t_stop(i))/(E_inc_stop(i).^3)).^2*(sigma_inc2(i))^2)+((1/(E_inc_stop(i).^2*2)).^2*Dt^2)+(t_stop(i)/(2*E_inc_stop(i)*Ang)).^2*dAng^2));
err_per_stop(i)=(sigma_dN_stop(i)/dN_dE_stop(i))*100;
end;

%% Energy spectrum
figure;
plot(E_inc_stop,dN_dE_stop/Ang)
xlabel('Released proton energy [MeV] ','FontSize',20);
ylabel('dN/dE \Omega [MeV^{-1}sr^{-1}]','FontSize',20);


% Number
NN_stop=trapz(E_inc_stop,dN_dE_stop)

